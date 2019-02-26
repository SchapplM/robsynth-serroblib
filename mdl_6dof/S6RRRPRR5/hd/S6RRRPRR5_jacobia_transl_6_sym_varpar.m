% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR5
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
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:17
% EndTime: 2019-02-26 22:18:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (209->43), mult. (195->55), div. (0->0), fcn. (210->10), ass. (0->38)
t35 = sin(qJ(5));
t67 = pkin(5) * t35 + qJ(4);
t34 = qJ(2) + qJ(3);
t31 = cos(t34);
t40 = -pkin(10) - pkin(9);
t64 = -pkin(3) - r_i_i_C(3);
t68 = (-t40 - t64) * t31;
t33 = qJ(5) + qJ(6);
t30 = cos(t33);
t58 = t30 * t31;
t28 = sin(t33);
t59 = t28 * t31;
t66 = r_i_i_C(1) * t59 + r_i_i_C(2) * t58 + t67 * t31;
t29 = sin(t34);
t32 = cos(qJ(2)) * pkin(2);
t65 = t67 * t29 + pkin(1) + t32 + t68;
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t54 = t39 * t30;
t5 = -t37 * t28 + t29 * t54;
t55 = t39 * t29;
t56 = t37 * t30;
t6 = t28 * t55 + t56;
t63 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t7 = t39 * t28 + t29 * t56;
t57 = t37 * t29;
t8 = -t28 * t57 + t54;
t62 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t38 = cos(qJ(5));
t60 = t38 * pkin(5);
t53 = pkin(4) + t60 + pkin(8) + pkin(7);
t47 = t64 * t29;
t45 = t66 * t37 + t40 * t57;
t44 = t66 * t39 + t40 * t55;
t43 = -sin(qJ(2)) * pkin(2) + t47;
t42 = t68 + (r_i_i_C(1) * t28 + r_i_i_C(2) * t30 + t67) * t29;
t17 = r_i_i_C(2) * t59;
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t65 * t37 + t53 * t39, t43 * t39 + t44, t39 * t47 + t44, t55 (-t35 * t37 + t38 * t55) * pkin(5) + t63, t63; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t53 * t37 + t65 * t39, t43 * t37 + t45, t37 * t47 + t45, t57 (t35 * t39 + t38 * t57) * pkin(5) + t62, t62; 0, t32 + t42, t42, -t31, t17 + (-r_i_i_C(1) * t30 - t60) * t31, -r_i_i_C(1) * t58 + t17;];
Ja_transl  = t1;
