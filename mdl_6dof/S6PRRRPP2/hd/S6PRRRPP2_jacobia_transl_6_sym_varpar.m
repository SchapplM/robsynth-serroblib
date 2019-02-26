% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_transl = S6PRRRPP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:09:32
% EndTime: 2019-02-26 20:09:33
% DurationCPUTime: 0.16s
% Computational Cost: add. (223->44), mult. (577->82), div. (0->0), fcn. (740->10), ass. (0->38)
t30 = sin(qJ(3));
t33 = cos(qJ(3));
t41 = pkin(9) - r_i_i_C(3) - qJ(6);
t52 = pkin(3) * t33 + t41 * t30 + pkin(2);
t28 = sin(pkin(6));
t51 = t28 * t30;
t50 = t28 * t33;
t34 = cos(qJ(2));
t49 = t28 * t34;
t29 = sin(qJ(4));
t48 = t29 * t33;
t32 = cos(qJ(4));
t47 = t32 * t33;
t46 = t33 * t34;
t45 = r_i_i_C(2) + qJ(5);
t44 = cos(pkin(6));
t43 = cos(pkin(10));
t42 = pkin(4) + pkin(5) + r_i_i_C(1);
t27 = sin(pkin(10));
t39 = t27 * t44;
t38 = t28 * t43;
t37 = t44 * t43;
t35 = t45 * t29 + t42 * t32 + pkin(3);
t31 = sin(qJ(2));
t24 = t44 * t30 + t31 * t50;
t23 = -t31 * t51 + t44 * t33;
t22 = -t31 * t39 + t43 * t34;
t21 = t43 * t31 + t34 * t39;
t20 = t27 * t34 + t31 * t37;
t19 = t27 * t31 - t34 * t37;
t13 = t24 * t29 + t32 * t49;
t12 = t22 * t33 + t27 * t51;
t11 = -t22 * t30 + t27 * t50;
t10 = t20 * t33 - t30 * t38;
t9 = -t20 * t30 - t33 * t38;
t3 = t12 * t29 - t21 * t32;
t1 = t10 * t29 - t19 * t32;
t2 = [0, t22 * pkin(8) + t45 * (-t21 * t48 - t22 * t32) + t42 * (-t21 * t47 + t22 * t29) - t52 * t21, t35 * t11 + t41 * t12, t45 * (t12 * t32 + t21 * t29) - t42 * t3, t3, t11; 0, t20 * pkin(8) + t45 * (-t19 * t48 - t20 * t32) + t42 * (-t19 * t47 + t20 * t29) - t52 * t19, t41 * t10 + t35 * t9, t45 * (t10 * t32 + t19 * t29) - t42 * t1, t1, t9; 1 (t45 * (t29 * t46 - t31 * t32) + t42 * (t29 * t31 + t32 * t46) + pkin(8) * t31 + t52 * t34) * t28, t35 * t23 + t41 * t24, t45 * (t24 * t32 - t29 * t49) - t42 * t13, t13, t23;];
Ja_transl  = t2;
