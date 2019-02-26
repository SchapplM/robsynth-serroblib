% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:46:21
% EndTime: 2019-02-26 22:46:22
% DurationCPUTime: 0.53s
% Computational Cost: add. (830->115), mult. (2302->198), div. (0->0), fcn. (3028->14), ass. (0->70)
t101 = cos(qJ(3));
t100 = sin(qJ(1));
t102 = cos(qJ(2));
t71 = sin(qJ(2));
t74 = cos(qJ(1));
t92 = cos(pkin(6));
t84 = t74 * t92;
t58 = t100 * t71 - t102 * t84;
t59 = t100 * t102 + t71 * t84;
t70 = sin(qJ(3));
t91 = cos(pkin(7));
t85 = t70 * t91;
t66 = sin(pkin(7));
t67 = sin(pkin(6));
t97 = t67 * t74;
t89 = t66 * t97;
t36 = t59 * t101 - t58 * t85 - t70 * t89;
t86 = t67 * t91;
t52 = t58 * t66 - t74 * t86;
t69 = sin(qJ(4));
t73 = cos(qJ(4));
t18 = t36 * t73 + t52 * t69;
t80 = t91 * t101;
t35 = t101 * t89 + t58 * t80 + t59 * t70;
t68 = sin(qJ(5));
t72 = cos(qJ(5));
t1 = t18 * t68 - t35 * t72;
t110 = t18 * t72 + t35 * t68;
t107 = -t36 * t69 + t52 * t73;
t81 = t92 * t100;
t77 = t102 * t81 + t74 * t71;
t87 = t67 * t100;
t106 = -t66 * t87 + t77 * t91;
t105 = r_i_i_C(1) + pkin(5);
t93 = r_i_i_C(3) + qJ(6);
t78 = t105 * t72 + t93 * t68 + pkin(4);
t104 = r_i_i_C(2) + pkin(12);
t103 = pkin(10) * t66;
t99 = t66 * t69;
t98 = t66 * t73;
t96 = t68 * t73;
t95 = t71 * t66;
t94 = t72 * t73;
t90 = t67 * t95;
t88 = t67 * t102;
t83 = t92 * t66;
t79 = -pkin(4) * t73 - t104 * t69 - pkin(3);
t75 = t100 * t86 + t77 * t66;
t60 = t74 * t102 - t71 * t81;
t57 = -t66 * t88 + t92 * t91;
t56 = (t102 * t101 - t71 * t85) * t67;
t55 = (t102 * t70 + t71 * t80) * t67;
t50 = t70 * t83 + (t101 * t71 + t102 * t85) * t67;
t49 = t67 * t71 * t70 - t101 * t83 - t80 * t88;
t46 = t56 * t73 + t69 * t90;
t44 = -t77 * t101 - t60 * t85;
t43 = t60 * t80 - t77 * t70;
t42 = -t58 * t101 - t59 * t85;
t41 = -t58 * t70 + t59 * t80;
t40 = t60 * t101 - t106 * t70;
t39 = t106 * t101 + t60 * t70;
t34 = t50 * t73 + t57 * t69;
t28 = t44 * t73 + t60 * t99;
t26 = t42 * t73 + t59 * t99;
t22 = t40 * t73 + t75 * t69;
t21 = t40 * t69 - t75 * t73;
t15 = t34 * t68 - t49 * t72;
t6 = t22 * t72 + t39 * t68;
t5 = t22 * t68 - t39 * t72;
t2 = [-t100 * pkin(1) - t59 * pkin(2) - t36 * pkin(3) - t18 * pkin(4) + pkin(9) * t97 - t52 * pkin(10) - t35 * pkin(11) - t93 * t1 + t104 * t107 - t105 * t110, t28 * pkin(4) + t44 * pkin(3) + t43 * pkin(11) - t77 * pkin(2) + t60 * t103 + t104 * (t44 * t69 - t60 * t98) + t105 * (t28 * t72 + t43 * t68) + t93 * (t28 * t68 - t43 * t72) t40 * pkin(11) + t93 * (-t39 * t96 - t40 * t72) + t105 * (-t39 * t94 + t40 * t68) + t79 * t39, t104 * t22 - t78 * t21, -t105 * t5 + t93 * t6, t5; t74 * pkin(1) + t60 * pkin(2) + t40 * pkin(3) + t22 * pkin(4) + pkin(9) * t87 + t75 * pkin(10) + t39 * pkin(11) + t104 * t21 + t105 * t6 + t93 * t5, t59 * t103 - t58 * pkin(2) + t42 * pkin(3) + t26 * pkin(4) + t41 * pkin(11) + t104 * (t42 * t69 - t59 * t98) + t105 * (t26 * t72 + t41 * t68) + t93 * (t26 * t68 - t41 * t72) t36 * pkin(11) + t105 * (-t35 * t94 + t36 * t68) + t93 * (-t35 * t96 - t36 * t72) + t79 * t35, t104 * t18 + t78 * t107, -t105 * t1 + t93 * t110, t1; 0, t56 * pkin(3) + t46 * pkin(4) + t55 * pkin(11) + (t102 * pkin(2) + pkin(10) * t95) * t67 + t104 * (t56 * t69 - t73 * t90) + t105 * (t46 * t72 + t55 * t68) + t93 * (t46 * t68 - t55 * t72) t50 * pkin(11) + t105 * (-t49 * t94 + t50 * t68) + t93 * (-t49 * t96 - t50 * t72) + t79 * t49, t104 * t34 + t78 * (-t50 * t69 + t57 * t73) t93 * (t34 * t72 + t49 * t68) - t105 * t15, t15;];
Ja_transl  = t2;
