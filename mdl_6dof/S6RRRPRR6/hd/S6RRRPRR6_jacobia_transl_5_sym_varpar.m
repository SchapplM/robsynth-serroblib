% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_transl = S6RRRPRR6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:47
% EndTime: 2019-02-26 22:18:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (157->32), mult. (138->45), div. (0->0), fcn. (153->10), ass. (0->27)
t21 = cos(qJ(2));
t19 = sin(qJ(2));
t30 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8);
t26 = t30 * t19;
t18 = qJ(3) + pkin(11);
t11 = pkin(4) * cos(t18) + cos(qJ(3)) * pkin(3);
t9 = pkin(2) + t11;
t34 = t21 * t9 + pkin(1) + t26;
t15 = qJ(5) + t18;
t13 = sin(t15);
t14 = cos(t15);
t22 = cos(qJ(1));
t20 = sin(qJ(1));
t29 = t20 * t21;
t5 = t13 * t29 + t14 * t22;
t6 = t13 * t22 - t14 * t29;
t33 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t28 = t21 * t22;
t7 = -t13 * t28 + t20 * t14;
t8 = t20 * t13 + t14 * t28;
t32 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t10 = pkin(4) * sin(t18) + sin(qJ(3)) * pkin(3);
t31 = pkin(7) + t10;
t25 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t24 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t23 = -t24 * t19 + t30 * t21;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t34 * t20 + t31 * t22, t23 * t22, -t10 * t28 + t20 * t11 + t32, t22 * t19, t32, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t31 * t20 + t34 * t22, t23 * t20, -t10 * t29 - t11 * t22 + t33, t20 * t19, t33, 0; 0, t24 * t21 + t26 (-t10 + t25) * t19, -t21, t25 * t19, 0;];
Ja_transl  = t1;
