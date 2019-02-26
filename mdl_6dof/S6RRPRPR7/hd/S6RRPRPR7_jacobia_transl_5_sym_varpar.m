% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:12
% EndTime: 2019-02-26 21:41:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (99->30), mult. (144->41), div. (0->0), fcn. (165->8), ass. (0->28)
t14 = sin(qJ(2));
t17 = cos(qJ(2));
t13 = sin(qJ(4));
t26 = pkin(4) * t13 + qJ(3);
t16 = cos(qJ(4));
t30 = pkin(4) * t16 + pkin(2) + pkin(3);
t20 = t26 * t14 + t30 * t17;
t31 = pkin(1) + t20;
t11 = qJ(4) + pkin(10);
t9 = sin(t11);
t29 = t17 * t9;
t18 = cos(qJ(1));
t28 = t18 * t14;
t27 = pkin(7) - r_i_i_C(3) - qJ(5) - pkin(8);
t15 = sin(qJ(1));
t10 = cos(t11);
t6 = t10 * t14 - t29;
t1 = t6 * t15;
t22 = t10 * t17 + t14 * t9;
t2 = t22 * t15;
t25 = t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
t3 = -t10 * t28 + t18 * t29;
t4 = t22 * t18;
t24 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
t23 = -t22 * r_i_i_C(1) - t6 * r_i_i_C(2);
t21 = pkin(4) * (-t13 * t17 + t14 * t16);
t19 = -t30 * t14 + t26 * t17;
t5 = [-t2 * r_i_i_C(1) - t1 * r_i_i_C(2) - t31 * t15 + t27 * t18, t19 * t18 - t24, t28, t18 * t21 + t24, -t15, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t27 * t15 + t31 * t18, t19 * t15 - t25, t15 * t14, t15 * t21 + t25, t18, 0; 0, t20 - t23, -t17 (-t13 * t14 - t16 * t17) * pkin(4) + t23, 0, 0;];
Ja_transl  = t5;
