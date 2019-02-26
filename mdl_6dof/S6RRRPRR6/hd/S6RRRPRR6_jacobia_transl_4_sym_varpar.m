% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function Ja_transl = S6RRRPRR6_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR6_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:47
% EndTime: 2019-02-26 22:18:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (77->29), mult. (106->44), div. (0->0), fcn. (117->8), ass. (0->24)
t14 = cos(qJ(2));
t11 = sin(qJ(2));
t24 = r_i_i_C(3) + qJ(4) + pkin(8);
t20 = t24 * t11;
t13 = cos(qJ(3));
t5 = pkin(3) * t13 + pkin(2);
t25 = t14 * t5 + pkin(1) + t20;
t10 = sin(qJ(3));
t23 = pkin(3) * t10;
t12 = sin(qJ(1));
t22 = t12 * t14;
t15 = cos(qJ(1));
t21 = t14 * t15;
t18 = pkin(7) + t23;
t8 = qJ(3) + pkin(11);
t6 = sin(t8);
t7 = cos(t8);
t17 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + t5;
t16 = -t17 * t11 + t24 * t14;
t4 = t12 * t6 + t7 * t21;
t3 = t12 * t7 - t6 * t21;
t2 = t15 * t6 - t7 * t22;
t1 = t15 * t7 + t6 * t22;
t9 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t25 * t12 + t18 * t15, t16 * t15, t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + (-t10 * t21 + t12 * t13) * pkin(3), t15 * t11, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t12 + t25 * t15, t16 * t12, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + (-t10 * t22 - t13 * t15) * pkin(3), t12 * t11, 0, 0; 0, t17 * t14 + t20 (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7 - t23) * t11, -t14, 0, 0;];
Ja_transl  = t9;
