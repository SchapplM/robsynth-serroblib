% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPP2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:09:22
% EndTime: 2019-02-26 20:09:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (174->44), mult. (453->82), div. (0->0), fcn. (579->10), ass. (0->35)
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t49 = pkin(9) + r_i_i_C(2);
t51 = pkin(3) * t32 + t49 * t29 + pkin(2);
t50 = pkin(4) + r_i_i_C(1);
t27 = sin(pkin(6));
t48 = t27 * t29;
t47 = t27 * t32;
t33 = cos(qJ(2));
t46 = t27 * t33;
t28 = sin(qJ(4));
t45 = t28 * t32;
t31 = cos(qJ(4));
t44 = t31 * t32;
t43 = t32 * t33;
t42 = r_i_i_C(3) + qJ(5);
t41 = cos(pkin(6));
t40 = cos(pkin(10));
t26 = sin(pkin(10));
t38 = t26 * t41;
t37 = t27 * t40;
t36 = t41 * t40;
t34 = t42 * t28 + t50 * t31 + pkin(3);
t30 = sin(qJ(2));
t24 = t41 * t29 + t30 * t47;
t22 = -t30 * t38 + t40 * t33;
t21 = t40 * t30 + t33 * t38;
t20 = t26 * t33 + t30 * t36;
t19 = t26 * t30 - t33 * t36;
t13 = t24 * t28 + t31 * t46;
t12 = t22 * t32 + t26 * t48;
t10 = t20 * t32 - t29 * t37;
t3 = t12 * t28 - t21 * t31;
t1 = t10 * t28 - t19 * t31;
t2 = [0, t22 * pkin(8) + t50 * (-t21 * t44 + t22 * t28) + t42 * (-t21 * t45 - t22 * t31) - t51 * t21, t49 * t12 + t34 * (-t22 * t29 + t26 * t47) t42 * (t12 * t31 + t21 * t28) - t50 * t3, t3, 0; 0, t20 * pkin(8) + t50 * (-t19 * t44 + t20 * t28) + t42 * (-t19 * t45 - t20 * t31) - t51 * t19, t49 * t10 + t34 * (-t20 * t29 - t32 * t37) t42 * (t10 * t31 + t19 * t28) - t50 * t1, t1, 0; 1 (t50 * (t28 * t30 + t31 * t43) + t42 * (t28 * t43 - t30 * t31) + pkin(8) * t30 + t51 * t33) * t27, t49 * t24 + t34 * (-t30 * t48 + t41 * t32) t42 * (t24 * t31 - t28 * t46) - t50 * t13, t13, 0;];
Ja_transl  = t2;
