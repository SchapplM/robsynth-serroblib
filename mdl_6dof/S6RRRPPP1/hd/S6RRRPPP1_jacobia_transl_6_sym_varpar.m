% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:47
% EndTime: 2019-02-26 22:02:47
% DurationCPUTime: 0.31s
% Computational Cost: add. (239->67), mult. (669->105), div. (0->0), fcn. (834->10), ass. (0->52)
t32 = sin(pkin(10));
t34 = cos(pkin(10));
t36 = sin(qJ(3));
t35 = cos(pkin(6));
t37 = sin(qJ(2));
t74 = t37 * t35;
t33 = sin(pkin(6));
t40 = cos(qJ(2));
t77 = t33 * t40;
t49 = t36 * t74 + t77;
t39 = cos(qJ(3));
t72 = t37 * t39;
t86 = -t32 * t49 + t34 * t72;
t57 = qJ(4) * t35 + pkin(9);
t54 = t57 * t37;
t84 = t40 * pkin(2) + pkin(1) + t54;
t41 = cos(qJ(1));
t66 = t41 * t39;
t38 = sin(qJ(1));
t70 = t38 * t36;
t22 = t40 * t70 + t66;
t67 = t41 * t36;
t69 = t39 * t40;
t23 = t38 * t69 - t67;
t73 = t37 * t38;
t62 = t33 * t73;
t83 = (-t22 * t35 + t62) * t34 - t23 * t32;
t44 = t32 * t72 + t49 * t34;
t64 = qJ(4) * t33;
t47 = pkin(3) * t39 + t36 * t64 + pkin(2);
t68 = t40 * t35;
t78 = t33 * t37;
t48 = t36 * t78 - t68;
t63 = pkin(4) + r_i_i_C(3) + qJ(6);
t65 = r_i_i_C(2) + qJ(5);
t81 = pkin(5) + r_i_i_C(1);
t82 = -t47 * t37 + t57 * t40 - t65 * t44 - t81 * t48 - t63 * t86;
t79 = t32 * t35;
t76 = t34 * t35;
t75 = t35 * t39;
t71 = t37 * t41;
t55 = (qJ(4) + t81) * t33;
t53 = t22 * t33 + t35 * t73;
t24 = t38 * t39 - t40 * t67;
t51 = t24 * t35 + t33 * t71;
t50 = t36 * t68 - t78;
t46 = t22 * t79 - t23 * t34 - t32 * t62;
t25 = t40 * t66 + t70;
t16 = -t24 * t33 + t35 * t71;
t4 = t25 * t34 + t51 * t32;
t3 = t25 * t32 - t51 * t34;
t1 = [-t23 * pkin(3) + t41 * pkin(8) - t22 * t64 - t84 * t38 + t63 * t46 - t81 * t53 + t65 * t83, t82 * t41, t24 * pkin(3) + t65 * (t24 * t32 + t25 * t76) + t63 * (t24 * t34 - t25 * t79) + t25 * t55, t16, t3, t4; t25 * pkin(3) + t38 * pkin(8) + t81 * t16 - t24 * t64 + t65 * t3 + t63 * t4 + t84 * t41, t82 * t38, -t22 * pkin(3) + t65 * (-t22 * t32 + t23 * t76) + t63 * (-t22 * t34 - t23 * t79) + t23 * t55, t53, -t83, -t46; 0, t54 + t81 * (t36 * t77 + t74) + t65 * (t32 * t69 + t50 * t34) + t47 * t40 + t63 * (-t32 * t50 + t34 * t69) (-t65 * (t32 * t36 - t34 * t75) - t63 * (t32 * t75 + t34 * t36) - pkin(3) * t36 + t39 * t55) * t37, t48, t44, t86;];
Ja_transl  = t1;
