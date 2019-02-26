% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:07
% EndTime: 2019-02-26 22:17:07
% DurationCPUTime: 0.21s
% Computational Cost: add. (303->41), mult. (364->59), div. (0->0), fcn. (436->10), ass. (0->37)
t32 = qJ(2) + qJ(3);
t29 = sin(t32);
t30 = cos(t32);
t52 = sin(qJ(5));
t53 = cos(qJ(5));
t16 = t29 * t52 + t30 * t53;
t35 = sin(qJ(1));
t10 = t16 * t35;
t57 = -r_i_i_C(3) - pkin(10);
t33 = sin(qJ(6));
t36 = cos(qJ(6));
t60 = t36 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(5);
t17 = t29 * t53 - t30 * t52;
t9 = t17 * t35;
t67 = -t57 * t10 + t60 * t9;
t37 = cos(qJ(1));
t46 = t37 * t52;
t47 = t37 * t53;
t11 = -t29 * t46 - t30 * t47;
t12 = -t29 * t47 + t30 * t46;
t66 = t57 * t11 - t60 * t12;
t65 = -t60 * t16 - t57 * t17;
t58 = pkin(3) + pkin(4);
t61 = t29 * qJ(4) + t58 * t30;
t31 = cos(qJ(2)) * pkin(2);
t59 = pkin(1) + t31 + t61;
t54 = -pkin(9) + pkin(8) + pkin(7);
t51 = qJ(4) * t30;
t50 = t58 * t29;
t44 = -t33 * r_i_i_C(1) - t36 * r_i_i_C(2);
t42 = -sin(qJ(2)) * pkin(2) - t50;
t41 = t35 * t51 - t67;
t40 = t37 * t51 - t66;
t39 = t61 - t65;
t2 = -t11 * t36 - t35 * t33;
t1 = t11 * t33 - t35 * t36;
t3 = [-t57 * t9 - t60 * t10 + (t44 + t54) * t37 - t59 * t35, t42 * t37 + t40, -t37 * t50 + t40, t37 * t29, t66, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -t11 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t57 * t12 + t54 * t35 + t59 * t37, t42 * t35 + t41, -t35 * t50 + t41, t35 * t29, t67 (-t10 * t33 + t37 * t36) * r_i_i_C(1) + (-t10 * t36 - t37 * t33) * r_i_i_C(2); 0, t31 + t39, t39, -t30, t65, t44 * t17;];
Ja_transl  = t3;
