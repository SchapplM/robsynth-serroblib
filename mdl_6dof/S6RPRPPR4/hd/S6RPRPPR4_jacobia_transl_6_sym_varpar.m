% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:53
% EndTime: 2019-02-26 20:40:53
% DurationCPUTime: 0.14s
% Computational Cost: add. (166->38), mult. (219->54), div. (0->0), fcn. (267->9), ass. (0->29)
t11 = pkin(9) + qJ(3);
t10 = cos(t11);
t24 = -r_i_i_C(3) - pkin(8) + qJ(4);
t9 = sin(t11);
t22 = t24 * t9;
t31 = t22 + t10 * pkin(3) + cos(pkin(9)) * pkin(2) + pkin(1);
t12 = sin(pkin(10));
t13 = cos(pkin(10));
t15 = sin(qJ(6));
t17 = cos(qJ(6));
t29 = pkin(4) + pkin(5);
t20 = t17 * r_i_i_C(1) - t15 * r_i_i_C(2) + t29;
t21 = t15 * r_i_i_C(1) + t17 * r_i_i_C(2) + qJ(5);
t30 = t21 * t12 + t20 * t13 + pkin(3);
t16 = sin(qJ(1));
t28 = t16 * t12;
t27 = t16 * t13;
t18 = cos(qJ(1));
t26 = t18 * t12;
t25 = t18 * t13;
t19 = t24 * t10 - t30 * t9;
t14 = -pkin(7) - qJ(2);
t6 = t10 * t25 + t28;
t5 = t10 * t26 - t27;
t4 = t10 * t27 - t26;
t3 = t10 * t28 + t25;
t2 = t5 * t15 + t6 * t17;
t1 = -t6 * t15 + t5 * t17;
t7 = [-t18 * t14 - t31 * t16 - t20 * t4 - t21 * t3, t16, t19 * t18, t18 * t9, t5, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(5) - t16 * t14 + t31 * t18 + t29 * t6, -t18, t19 * t16, t16 * t9, t3 (-t4 * t15 + t3 * t17) * r_i_i_C(1) + (-t3 * t15 - t4 * t17) * r_i_i_C(2); 0, 0, t30 * t10 + t22, -t10, t9 * t12 ((t12 * t17 - t13 * t15) * r_i_i_C(1) + (-t12 * t15 - t13 * t17) * r_i_i_C(2)) * t9;];
Ja_transl  = t7;
