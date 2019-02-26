% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR6
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
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:47
% EndTime: 2019-02-26 22:18:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (284->40), mult. (182->54), div. (0->0), fcn. (201->12), ass. (0->33)
t25 = cos(qJ(2));
t23 = sin(qJ(2));
t36 = r_i_i_C(3) + pkin(10) + pkin(9) + qJ(4) + pkin(8);
t30 = t36 * t23;
t22 = qJ(3) + pkin(11);
t19 = qJ(5) + t22;
t17 = cos(t19);
t11 = pkin(5) * t17 + pkin(4) * cos(t22) + cos(qJ(3)) * pkin(3);
t9 = pkin(2) + t11;
t41 = t25 * t9 + pkin(1) + t30;
t18 = qJ(6) + t19;
t13 = sin(t18);
t14 = cos(t18);
t26 = cos(qJ(1));
t32 = t26 * t14;
t24 = sin(qJ(1));
t35 = t24 * t25;
t5 = t13 * t35 + t32;
t33 = t26 * t13;
t6 = -t14 * t35 + t33;
t40 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = t24 * t14 - t25 * t33;
t8 = t24 * t13 + t25 * t32;
t39 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t16 = sin(t19);
t38 = pkin(5) * t16;
t10 = t38 + pkin(4) * sin(t22) + sin(qJ(3)) * pkin(3);
t37 = pkin(7) + t10;
t34 = t25 * t26;
t29 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t28 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t27 = -t28 * t23 + t36 * t25;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t41 * t24 + t37 * t26, t27 * t26, -t10 * t34 + t24 * t11 + t39, t26 * t23 (-t16 * t34 + t17 * t24) * pkin(5) + t39, t39; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t37 * t24 + t41 * t26, t27 * t24, -t10 * t35 - t26 * t11 + t40, t24 * t23 (-t16 * t35 - t17 * t26) * pkin(5) + t40, t40; 0, t28 * t25 + t30 (-t10 + t29) * t23, -t25 (t29 - t38) * t23, t29 * t23;];
Ja_transl  = t1;
