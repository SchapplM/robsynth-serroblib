% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR12
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR12_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR12_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobia_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:07
% EndTime: 2019-02-26 21:21:08
% DurationCPUTime: 0.30s
% Computational Cost: add. (260->63), mult. (734->112), div. (0->0), fcn. (961->14), ass. (0->49)
t29 = sin(qJ(4));
t32 = cos(qJ(4));
t28 = cos(pkin(6));
t25 = cos(pkin(14));
t34 = cos(qJ(1));
t46 = t34 * t25;
t21 = sin(pkin(14));
t31 = sin(qJ(1));
t51 = t31 * t21;
t15 = -t28 * t46 + t51;
t23 = sin(pkin(7));
t27 = cos(pkin(7));
t24 = sin(pkin(6));
t47 = t34 * t24;
t11 = -t15 * t23 + t27 * t47;
t22 = sin(pkin(8));
t26 = cos(pkin(8));
t48 = t34 * t21;
t49 = t31 * t25;
t16 = t28 * t48 + t49;
t30 = sin(qJ(3));
t33 = cos(qJ(3));
t44 = t23 * t47;
t5 = (t15 * t27 + t44) * t33 + t16 * t30;
t41 = t11 * t22 + t26 * t5;
t52 = t27 * t30;
t6 = t15 * t52 - t16 * t33 + t30 * t44;
t62 = t29 * t41 + t6 * t32;
t61 = -t6 * t29 + t32 * t41;
t58 = pkin(11) + r_i_i_C(3);
t55 = t23 * t28;
t54 = t26 * t29;
t53 = t26 * t32;
t50 = t31 * t24;
t45 = t24 * qJ(2);
t43 = t58 * t22;
t17 = -t28 * t49 - t48;
t13 = -t17 * t23 + t27 * t50;
t18 = -t28 * t51 + t46;
t35 = t17 * t27 + t23 * t50;
t7 = -t18 * t30 + t33 * t35;
t38 = t13 * t22 + t26 * t7;
t9 = t33 * t55 + (t25 * t27 * t33 - t21 * t30) * t24;
t37 = (-t23 * t24 * t25 + t27 * t28) * t22 + t26 * t9;
t10 = t30 * t55 + (t21 * t33 + t25 * t52) * t24;
t8 = t18 * t33 + t30 * t35;
t2 = t29 * t38 + t8 * t32;
t1 = -t8 * t29 + t32 * t38;
t3 = [t62 * r_i_i_C(1) + t61 * r_i_i_C(2) + t6 * pkin(3) - t16 * pkin(2) - t31 * pkin(1) + t34 * t45 + t11 * pkin(10) + t58 * (t11 * t26 - t22 * t5) t50 (t32 * t7 - t54 * t8) * r_i_i_C(1) + (-t29 * t7 - t53 * t8) * r_i_i_C(2) + t7 * pkin(3) + t8 * t43, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; t34 * pkin(1) + t18 * pkin(2) + t8 * pkin(3) + t13 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t31 * t45 + t58 * (t13 * t26 - t22 * t7) -t47 (-t32 * t5 + t54 * t6) * r_i_i_C(1) + (t29 * t5 + t53 * t6) * r_i_i_C(2) - t5 * pkin(3) - t6 * t43, -r_i_i_C(1) * t61 + r_i_i_C(2) * t62, 0, 0; 0, t28 (r_i_i_C(1) * t32 - r_i_i_C(2) * t29 + pkin(3)) * t9 + ((-r_i_i_C(1) * t29 - r_i_i_C(2) * t32) * t26 + t43) * t10 (-t10 * t29 + t32 * t37) * r_i_i_C(1) + (-t10 * t32 - t29 * t37) * r_i_i_C(2), 0, 0;];
Ja_transl  = t3;
