% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP3
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
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPP3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:01
% EndTime: 2019-02-26 20:10:01
% DurationCPUTime: 0.19s
% Computational Cost: add. (174->44), mult. (453->82), div. (0->0), fcn. (579->10), ass. (0->35)
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t47 = pkin(9) + r_i_i_C(1);
t49 = pkin(3) * t30 + t47 * t27 + pkin(2);
t48 = pkin(4) - r_i_i_C(2);
t25 = sin(pkin(6));
t46 = t25 * t27;
t45 = t25 * t30;
t31 = cos(qJ(2));
t44 = t25 * t31;
t26 = sin(qJ(4));
t43 = t26 * t30;
t29 = cos(qJ(4));
t42 = t29 * t30;
t41 = t30 * t31;
t40 = r_i_i_C(3) + qJ(5);
t39 = cos(pkin(6));
t38 = cos(pkin(10));
t24 = sin(pkin(10));
t36 = t24 * t39;
t35 = t25 * t38;
t34 = t39 * t38;
t32 = t40 * t26 + t48 * t29 + pkin(3);
t28 = sin(qJ(2));
t22 = t39 * t27 + t28 * t45;
t20 = -t28 * t36 + t38 * t31;
t19 = t38 * t28 + t31 * t36;
t18 = t24 * t31 + t28 * t34;
t17 = t24 * t28 - t31 * t34;
t13 = t22 * t26 + t29 * t44;
t12 = t20 * t30 + t24 * t46;
t10 = t18 * t30 - t27 * t35;
t3 = t12 * t26 - t19 * t29;
t1 = t10 * t26 - t17 * t29;
t2 = [0, t20 * pkin(8) + t48 * (-t19 * t42 + t20 * t26) + t40 * (-t19 * t43 - t20 * t29) - t49 * t19, t47 * t12 + t32 * (-t20 * t27 + t24 * t45) t40 * (t12 * t29 + t19 * t26) - t48 * t3, t3, 0; 0, t18 * pkin(8) + t48 * (-t17 * t42 + t18 * t26) + t40 * (-t17 * t43 - t18 * t29) - t49 * t17, t47 * t10 + t32 * (-t18 * t27 - t30 * t35) t40 * (t10 * t29 + t17 * t26) - t48 * t1, t1, 0; 1 (t48 * (t26 * t28 + t29 * t41) + t40 * (t26 * t41 - t28 * t29) + pkin(8) * t28 + t49 * t31) * t25, t47 * t22 + t32 * (-t28 * t46 + t39 * t30) t40 * (t22 * t29 - t26 * t44) - t48 * t13, t13, 0;];
Ja_transl  = t2;
