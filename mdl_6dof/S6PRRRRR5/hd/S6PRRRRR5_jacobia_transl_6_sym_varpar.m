% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:10
% EndTime: 2019-02-26 20:21:10
% DurationCPUTime: 0.31s
% Computational Cost: add. (597->81), mult. (1489->149), div. (0->0), fcn. (1942->16), ass. (0->59)
t46 = sin(pkin(7));
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t48 = cos(pkin(13));
t75 = cos(pkin(6));
t69 = t48 * t75;
t73 = sin(pkin(13));
t60 = t73 * t52 - t55 * t69;
t74 = cos(pkin(7));
t47 = sin(pkin(6));
t78 = t47 * t48;
t88 = t46 * t78 + t60 * t74;
t64 = t75 * t73;
t61 = t48 * t52 + t55 * t64;
t70 = t47 * t73;
t87 = -t46 * t70 + t61 * t74;
t36 = t52 * t69 + t73 * t55;
t51 = sin(qJ(3));
t81 = cos(qJ(3));
t15 = t36 * t51 + t88 * t81;
t45 = qJ(5) + qJ(6);
t43 = sin(t45);
t44 = cos(t45);
t16 = t36 * t81 - t88 * t51;
t29 = t60 * t46 - t74 * t78;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t8 = t16 * t54 + t29 * t50;
t86 = (t15 * t44 - t8 * t43) * r_i_i_C(1) + (-t15 * t43 - t8 * t44) * r_i_i_C(2);
t37 = t48 * t55 - t52 * t64;
t18 = t37 * t81 - t87 * t51;
t30 = t61 * t46 + t74 * t70;
t10 = t18 * t54 + t30 * t50;
t17 = t37 * t51 + t87 * t81;
t85 = (-t10 * t43 + t17 * t44) * r_i_i_C(1) + (-t10 * t44 - t17 * t43) * r_i_i_C(2);
t67 = t75 * t46;
t68 = t51 * t74;
t28 = t51 * t67 + (t81 * t52 + t55 * t68) * t47;
t76 = t47 * t55;
t35 = -t46 * t76 + t75 * t74;
t20 = t28 * t54 + t35 * t50;
t65 = t74 * t81;
t77 = t47 * t52;
t27 = t51 * t77 - t65 * t76 - t81 * t67;
t84 = (-t20 * t43 + t27 * t44) * r_i_i_C(1) + (-t20 * t44 - t27 * t43) * r_i_i_C(2);
t83 = pkin(9) * t46;
t82 = r_i_i_C(3) + pkin(12) + pkin(11);
t80 = t46 * t50;
t79 = t46 * t54;
t72 = t46 * t77;
t53 = cos(qJ(5));
t63 = t53 * pkin(5) + r_i_i_C(1) * t44 - r_i_i_C(2) * t43 + pkin(4);
t49 = sin(qJ(5));
t62 = t49 * pkin(5) + t43 * r_i_i_C(1) + t44 * r_i_i_C(2) + pkin(10);
t59 = -t82 * t50 - t63 * t54 - pkin(3);
t34 = (-t52 * t68 + t81 * t55) * t47;
t24 = -t37 * t68 - t61 * t81;
t22 = -t36 * t68 - t60 * t81;
t1 = [0, t24 * pkin(3) - t61 * pkin(2) + t37 * t83 + t63 * (t24 * t54 + t37 * t80) + t62 * (t37 * t65 - t61 * t51) + t82 * (t24 * t50 - t37 * t79) t59 * t17 + t62 * t18, t82 * t10 + t63 * (-t18 * t50 + t30 * t54) (-t10 * t49 + t17 * t53) * pkin(5) + t85, t85; 0, t22 * pkin(3) - t60 * pkin(2) + t36 * t83 + t63 * (t22 * t54 + t36 * t80) + t62 * (t36 * t65 - t60 * t51) + t82 * (t22 * t50 - t36 * t79) t59 * t15 + t62 * t16, t82 * t8 + t63 * (-t16 * t50 + t29 * t54) (t15 * t53 - t49 * t8) * pkin(5) + t86, t86; 1, t34 * pkin(3) + t63 * (t34 * t54 + t50 * t72) + t82 * (t34 * t50 - t54 * t72) + (t55 * pkin(2) + t52 * t83 + t62 * (t51 * t55 + t52 * t65)) * t47, t59 * t27 + t62 * t28, t82 * t20 + t63 * (-t28 * t50 + t35 * t54) (-t20 * t49 + t27 * t53) * pkin(5) + t84, t84;];
Ja_transl  = t1;
