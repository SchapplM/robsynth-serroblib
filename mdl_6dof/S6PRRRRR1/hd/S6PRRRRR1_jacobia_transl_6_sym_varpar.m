% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:36
% EndTime: 2019-02-26 20:18:36
% DurationCPUTime: 0.18s
% Computational Cost: add. (459->51), mult. (515->80), div. (0->0), fcn. (638->14), ass. (0->39)
t61 = pkin(11) + r_i_i_C(3);
t37 = sin(qJ(6));
t39 = cos(qJ(6));
t60 = t39 * r_i_i_C(1) - t37 * r_i_i_C(2) + pkin(5);
t33 = qJ(3) + qJ(4);
t29 = cos(t33);
t23 = pkin(4) * t29 + cos(qJ(3)) * pkin(3);
t30 = qJ(5) + t33;
t26 = sin(t30);
t27 = cos(t30);
t59 = t61 * t26 + t60 * t27 + pkin(2) + t23;
t34 = sin(pkin(12));
t35 = sin(pkin(6));
t55 = t34 * t35;
t38 = sin(qJ(2));
t54 = t34 * t38;
t40 = cos(qJ(2));
t53 = t34 * t40;
t52 = t35 * t38;
t51 = t35 * t40;
t50 = cos(pkin(12));
t49 = t35 * t50;
t48 = t50 * t38;
t47 = t50 * t40;
t45 = t37 * r_i_i_C(1) + t39 * r_i_i_C(2) + pkin(8) + pkin(9) + pkin(10);
t36 = cos(pkin(6));
t17 = t36 * t48 + t53;
t8 = t17 * t27 - t26 * t49;
t44 = t61 * t8 + t60 * (-t17 * t26 - t27 * t49);
t19 = -t36 * t54 + t47;
t10 = t19 * t27 + t26 * t55;
t43 = t61 * t10 + t60 * (-t19 * t26 + t27 * t55);
t15 = t36 * t26 + t27 * t52;
t42 = t61 * t15 + t60 * (-t26 * t52 + t36 * t27);
t28 = sin(t33);
t22 = -pkin(4) * t28 - sin(qJ(3)) * pkin(3);
t18 = t36 * t53 + t48;
t16 = -t36 * t47 + t54;
t1 = [0, -t18 * t59 + t45 * t19, t19 * t22 + t23 * t55 + t43 (-t19 * t28 + t29 * t55) * pkin(4) + t43, t43 (-t10 * t37 + t18 * t39) * r_i_i_C(1) + (-t10 * t39 - t18 * t37) * r_i_i_C(2); 0, -t16 * t59 + t45 * t17, t17 * t22 - t23 * t49 + t44 (-t17 * t28 - t29 * t49) * pkin(4) + t44, t44 (t16 * t39 - t8 * t37) * r_i_i_C(1) + (-t16 * t37 - t8 * t39) * r_i_i_C(2); 1 (t45 * t38 + t59 * t40) * t35, t22 * t52 + t36 * t23 + t42 (-t28 * t52 + t29 * t36) * pkin(4) + t42, t42 (-t15 * t37 - t39 * t51) * r_i_i_C(1) + (-t15 * t39 + t37 * t51) * r_i_i_C(2);];
Ja_transl  = t1;
