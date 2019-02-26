% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:00
% EndTime: 2019-02-26 22:45:00
% DurationCPUTime: 0.26s
% Computational Cost: add. (435->67), mult. (802->110), div. (0->0), fcn. (1027->12), ass. (0->47)
t78 = pkin(5) + r_i_i_C(1);
t67 = r_i_i_C(3) + qJ(6);
t49 = sin(qJ(2));
t50 = sin(qJ(1));
t53 = cos(qJ(2));
t66 = cos(pkin(6));
t76 = cos(qJ(1));
t60 = t66 * t76;
t35 = t49 * t60 + t50 * t53;
t48 = sin(qJ(3));
t52 = cos(qJ(3));
t46 = sin(pkin(6));
t64 = t46 * t76;
t23 = t35 * t52 - t48 * t64;
t34 = t50 * t49 - t53 * t60;
t45 = qJ(4) + qJ(5);
t43 = sin(t45);
t44 = cos(t45);
t7 = t23 * t43 - t34 * t44;
t80 = t23 * t44 + t34 * t43;
t51 = cos(qJ(4));
t42 = t51 * pkin(4) + pkin(3);
t77 = r_i_i_C(2) + pkin(11) + pkin(10);
t79 = t42 * t52 + t77 * t48 + pkin(2);
t55 = t67 * t43 + t78 * t44 + t42;
t73 = t43 * t52;
t72 = t44 * t52;
t71 = t46 * t50;
t70 = t46 * t52;
t69 = t46 * t53;
t68 = t52 * t53;
t47 = sin(qJ(4));
t65 = pkin(4) * t47 + pkin(9);
t62 = t50 * t66;
t61 = t67 * t80 - t78 * t7;
t37 = -t49 * t62 + t76 * t53;
t27 = t37 * t52 + t48 * t71;
t36 = t76 * t49 + t53 * t62;
t11 = t27 * t43 - t36 * t44;
t12 = t27 * t44 + t36 * t43;
t59 = -t78 * t11 + t67 * t12;
t33 = t66 * t48 + t49 * t70;
t20 = t33 * t43 + t44 * t69;
t58 = t67 * (t33 * t44 - t43 * t69) - t78 * t20;
t56 = -t35 * t48 - t52 * t64;
t26 = t37 * t48 - t50 * t70;
t1 = [-t50 * pkin(1) - t35 * pkin(2) + pkin(8) * t64 - t23 * t42 - t65 * t34 + t77 * t56 - t67 * t7 - t78 * t80, t65 * t37 + t78 * (-t36 * t72 + t37 * t43) + t67 * (-t36 * t73 - t37 * t44) - t79 * t36, -t55 * t26 + t77 * t27 (-t27 * t47 + t36 * t51) * pkin(4) + t59, t59, t11; t76 * pkin(1) + t37 * pkin(2) + pkin(8) * t71 + t67 * t11 + t78 * t12 + t77 * t26 + t27 * t42 + t65 * t36, t65 * t35 + t78 * (-t34 * t72 + t35 * t43) + t67 * (-t34 * t73 - t35 * t44) - t79 * t34, t77 * t23 + t55 * t56 (-t23 * t47 + t34 * t51) * pkin(4) + t61, t61, t7; 0 (t78 * (t43 * t49 + t44 * t68) + t67 * (t43 * t68 - t44 * t49) + t65 * t49 + t79 * t53) * t46, t77 * t33 + t55 * (-t46 * t49 * t48 + t66 * t52) (-t33 * t47 - t51 * t69) * pkin(4) + t58, t58, t20;];
Ja_transl  = t1;
