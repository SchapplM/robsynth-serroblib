% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR9_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR9_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobia_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:27
% EndTime: 2019-02-26 22:20:27
% DurationCPUTime: 0.18s
% Computational Cost: add. (172->58), mult. (443->100), div. (0->0), fcn. (561->12), ass. (0->39)
t24 = sin(qJ(3));
t43 = t24 * pkin(3);
t27 = cos(qJ(3));
t42 = t27 * pkin(3);
t21 = sin(pkin(6));
t26 = sin(qJ(1));
t41 = t21 * t26;
t28 = cos(qJ(2));
t40 = t21 * t28;
t29 = cos(qJ(1));
t39 = t21 * t29;
t38 = pkin(10) + qJ(4);
t37 = cos(pkin(6));
t36 = t26 * t37;
t35 = t29 * t37;
t20 = sin(pkin(7));
t23 = cos(pkin(7));
t19 = sin(pkin(13));
t22 = cos(pkin(13));
t32 = t27 * t19 + t24 * t22;
t4 = t32 * t20;
t34 = r_i_i_C(1) * t4 + t20 * t43 + t38 * t23 + pkin(9);
t25 = sin(qJ(2));
t11 = -t29 * t25 - t28 * t36;
t12 = -t25 * t36 + t29 * t28;
t13 = t24 * t19 - t27 * t22;
t6 = t32 * t23;
t33 = t11 * t6 - t12 * t13;
t18 = pkin(2) + t42;
t31 = t13 * r_i_i_C(1) + r_i_i_C(2) * t32 - t18;
t5 = t13 * t23;
t8 = -t38 * t20 + t23 * t43;
t30 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2) - t20 * r_i_i_C(3) + t8;
t10 = t25 * t35 + t26 * t28;
t9 = t26 * t25 - t28 * t35;
t3 = t13 * t20;
t2 = -t11 * t20 + t23 * t41;
t1 = -t11 * t5 - t12 * t32 - t3 * t41;
t7 = [-t26 * pkin(1) + t31 * t10 + t30 * t9 + (-r_i_i_C(2) * t3 + r_i_i_C(3) * t23 + t34) * t39, -t31 * t11 - t30 * t12, t1 * r_i_i_C(1) + (-t4 * t41 - t33) * r_i_i_C(2) + (-t12 * t24 + (t11 * t23 + t20 * t41) * t27) * pkin(3), t2, 0, 0; t29 * pkin(1) + t33 * r_i_i_C(1) + t1 * r_i_i_C(2) + t2 * r_i_i_C(3) + t11 * t8 + t12 * t18 + t34 * t41, -t30 * t10 + t31 * t9 (-t10 * t32 + t3 * t39 + t9 * t5) * r_i_i_C(1) + (t10 * t13 + t4 * t39 + t9 * t6) * r_i_i_C(2) + (-t10 * t24 + (-t20 * t39 - t23 * t9) * t27) * pkin(3), t9 * t20 - t23 * t39, 0, 0; 0 (-t30 * t25 - t31 * t28) * t21 (-t37 * t3 + (-t25 * t32 - t28 * t5) * t21) * r_i_i_C(1) + (-t37 * t4 + (t13 * t25 - t28 * t6) * t21) * r_i_i_C(2) - t21 * t25 * t43 + (t37 * t20 + t23 * t40) * t42, -t20 * t40 + t37 * t23, 0, 0;];
Ja_transl  = t7;
