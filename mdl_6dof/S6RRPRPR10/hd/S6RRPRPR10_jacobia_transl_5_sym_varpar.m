% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR10_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR10_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:00
% EndTime: 2019-02-26 21:43:00
% DurationCPUTime: 0.13s
% Computational Cost: add. (177->35), mult. (273->56), div. (0->0), fcn. (342->10), ass. (0->30)
t15 = cos(pkin(11)) * pkin(3) + pkin(2);
t18 = pkin(11) + qJ(4);
t16 = sin(t18);
t17 = cos(t18);
t32 = r_i_i_C(3) + qJ(5);
t37 = pkin(4) - r_i_i_C(2);
t38 = t32 * t16 + t37 * t17 + t15;
t36 = r_i_i_C(1) + pkin(9) + qJ(3);
t20 = sin(pkin(6));
t22 = sin(qJ(2));
t35 = t20 * t22;
t23 = sin(qJ(1));
t34 = t20 * t23;
t25 = cos(qJ(1));
t33 = t20 * t25;
t31 = cos(pkin(6));
t30 = t23 * t31;
t29 = t25 * t31;
t28 = t20 * (pkin(3) * sin(pkin(11)) + pkin(8));
t24 = cos(qJ(2));
t10 = t22 * t29 + t23 * t24;
t1 = t10 * t16 + t17 * t33;
t27 = -t10 * t17 + t16 * t33;
t12 = -t22 * t30 + t25 * t24;
t11 = t25 * t22 + t24 * t30;
t9 = t23 * t22 - t24 * t29;
t7 = t16 * t35 - t31 * t17;
t6 = t12 * t17 + t16 * t34;
t5 = t12 * t16 - t17 * t34;
t2 = [-t23 * pkin(1) - t32 * t1 - t10 * t15 + t25 * t28 + t37 * t27 - t36 * t9, -t11 * t38 + t36 * t12, t11, t32 * t6 - t37 * t5, t5, 0; t25 * pkin(1) + t36 * t11 + t12 * t15 + t23 * t28 + t32 * t5 + t37 * t6, t36 * t10 - t38 * t9, t9, -t37 * t1 - t32 * t27, t1, 0; 0 (t36 * t22 + t38 * t24) * t20, -t20 * t24, t32 * (t31 * t16 + t17 * t35) - t37 * t7, t7, 0;];
Ja_transl  = t2;
