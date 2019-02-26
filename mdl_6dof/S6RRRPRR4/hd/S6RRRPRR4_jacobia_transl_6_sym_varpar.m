% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:42
% EndTime: 2019-02-26 22:17:42
% DurationCPUTime: 0.16s
% Computational Cost: add. (256->43), mult. (177->56), div. (0->0), fcn. (192->12), ass. (0->40)
t28 = pkin(11) + qJ(5);
t22 = cos(t28);
t12 = pkin(5) * t22 + cos(pkin(11)) * pkin(4) + pkin(3);
t29 = qJ(2) + qJ(3);
t24 = sin(t29);
t25 = cos(t29);
t27 = -pkin(10) - pkin(9) - qJ(4);
t56 = t25 * t12 + (r_i_i_C(3) - t27) * t24;
t26 = cos(qJ(2)) * pkin(2);
t55 = pkin(1) + t26 + t56;
t23 = qJ(6) + t28;
t18 = cos(t23);
t32 = cos(qJ(1));
t43 = t32 * t18;
t17 = sin(t23);
t31 = sin(qJ(1));
t46 = t31 * t17;
t5 = t25 * t46 + t43;
t44 = t32 * t17;
t45 = t31 * t18;
t6 = -t25 * t45 + t44;
t54 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t25 * t44 + t45;
t8 = t25 * t43 + t46;
t53 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t21 = sin(t28);
t52 = pkin(5) * t21;
t51 = r_i_i_C(1) * t18;
t50 = r_i_i_C(2) * t17;
t48 = t25 * t31;
t47 = t25 * t32;
t39 = t24 * t50;
t42 = r_i_i_C(3) * t48 + t31 * t39;
t41 = r_i_i_C(3) * t47 + t32 * t39;
t40 = t52 + sin(pkin(11)) * pkin(4) + pkin(8) + pkin(7);
t37 = -r_i_i_C(1) * t17 - r_i_i_C(2) * t18;
t36 = -t25 * t27 + (-t12 - t51) * t24;
t35 = (-t50 + t51) * t25 + t56;
t34 = -sin(qJ(2)) * pkin(2) + t36;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t55 * t31 + t40 * t32, t32 * t34 + t41, t32 * t36 + t41, t32 * t24 (-t21 * t47 + t22 * t31) * pkin(5) + t53, t53; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t40 * t31 + t55 * t32, t31 * t34 + t42, t31 * t36 + t42, t31 * t24 (-t21 * t48 - t22 * t32) * pkin(5) + t54, t54; 0, t26 + t35, t35, -t25 (t37 - t52) * t24, t37 * t24;];
Ja_transl  = t1;
