% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:42
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.23s
% Computational Cost: add. (498->55), mult. (1136->97), div. (0->0), fcn. (1492->16), ass. (0->52)
t76 = pkin(11) + r_i_i_C(3);
t38 = sin(qJ(6));
t41 = cos(qJ(6));
t75 = t41 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
t37 = cos(pkin(12));
t63 = sin(pkin(13));
t64 = sin(pkin(12));
t54 = t64 * t63;
t66 = cos(pkin(13));
t61 = t37 * t66;
t68 = cos(pkin(6));
t46 = -t68 * t61 + t54;
t36 = sin(pkin(6));
t65 = sin(pkin(7));
t62 = t36 * t65;
t67 = cos(pkin(7));
t74 = t37 * t62 + t46 * t67;
t55 = t64 * t66;
t60 = t37 * t63;
t47 = t68 * t55 + t60;
t59 = t64 * t36;
t73 = t47 * t67 - t65 * t59;
t70 = cos(qJ(3));
t69 = t37 * t36;
t57 = t68 * t65;
t56 = t67 * t66;
t53 = t38 * r_i_i_C(1) + t41 * r_i_i_C(2) + pkin(9) + pkin(10);
t27 = t68 * t60 + t55;
t40 = sin(qJ(3));
t17 = t27 * t70 - t74 * t40;
t22 = t46 * t65 - t67 * t69;
t35 = qJ(4) + qJ(5);
t33 = sin(t35);
t34 = cos(t35);
t8 = t17 * t34 + t22 * t33;
t51 = t76 * t8 + t75 * (-t17 * t33 + t22 * t34);
t28 = -t68 * t54 + t61;
t19 = t28 * t70 - t73 * t40;
t23 = t47 * t65 + t67 * t59;
t10 = t19 * t34 + t23 * t33;
t50 = t76 * t10 + t75 * (-t19 * t33 + t23 * t34);
t21 = t40 * t57 + (t40 * t56 + t70 * t63) * t36;
t26 = -t66 * t62 + t68 * t67;
t15 = t21 * t34 + t26 * t33;
t49 = t76 * t15 + t75 * (-t21 * t33 + t26 * t34);
t42 = cos(qJ(4));
t48 = -t42 * pkin(4) - t76 * t33 - t75 * t34 - pkin(3);
t39 = sin(qJ(4));
t20 = -t70 * t57 + (t40 * t63 - t56 * t70) * t36;
t18 = t28 * t40 + t73 * t70;
t16 = t27 * t40 + t74 * t70;
t1 = [0, t59, t48 * t18 + t53 * t19 (-t19 * t39 + t23 * t42) * pkin(4) + t50, t50 (-t10 * t38 + t18 * t41) * r_i_i_C(1) + (-t10 * t41 - t18 * t38) * r_i_i_C(2); 0, -t69, t48 * t16 + t53 * t17 (-t17 * t39 + t22 * t42) * pkin(4) + t51, t51 (t16 * t41 - t8 * t38) * r_i_i_C(1) + (-t16 * t38 - t8 * t41) * r_i_i_C(2); 1, t68, t48 * t20 + t53 * t21 (-t21 * t39 + t26 * t42) * pkin(4) + t49, t49 (-t15 * t38 + t20 * t41) * r_i_i_C(1) + (-t15 * t41 - t20 * t38) * r_i_i_C(2);];
Ja_transl  = t1;
