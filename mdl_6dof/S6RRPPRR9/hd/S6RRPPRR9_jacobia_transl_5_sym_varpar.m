% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:12
% EndTime: 2019-02-26 21:33:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (88->32), mult. (208->54), div. (0->0), fcn. (260->8), ass. (0->26)
t10 = sin(pkin(6));
t12 = sin(qJ(2));
t27 = t10 * t12;
t13 = sin(qJ(1));
t26 = t10 * t13;
t14 = cos(qJ(5));
t25 = t10 * t14;
t16 = cos(qJ(1));
t24 = t10 * t16;
t23 = pkin(2) + qJ(4);
t22 = cos(pkin(6));
t21 = pkin(3) + pkin(4) + pkin(8);
t20 = r_i_i_C(3) + pkin(9) - qJ(3);
t19 = t13 * t22;
t18 = t16 * t22;
t11 = sin(qJ(5));
t17 = t11 * r_i_i_C(1) + t14 * r_i_i_C(2) + t23;
t15 = cos(qJ(2));
t8 = t14 * t24;
t6 = -t12 * t19 + t15 * t16;
t5 = t16 * t12 + t15 * t19;
t4 = t12 * t18 + t13 * t15;
t3 = t12 * t13 - t15 * t18;
t2 = t11 * t6 + t13 * t25;
t1 = -t11 * t26 + t14 * t6;
t7 = [-t13 * pkin(1) + t8 * r_i_i_C(1) - t17 * t4 + t20 * t3 + (-t11 * r_i_i_C(2) + t21) * t24, -t17 * t5 - t20 * t6, t5, t6, r_i_i_C(1) * t1 - t2 * r_i_i_C(2), 0; pkin(1) * t16 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t20 * t5 + t21 * t26 + t23 * t6, -t17 * t3 - t20 * t4, t3, t4 (t11 * t24 + t4 * t14) * r_i_i_C(1) + (-t11 * t4 + t8) * r_i_i_C(2), 0; 0 (-t20 * t12 + t17 * t15) * t10, -t10 * t15, t27 (-t22 * t11 + t12 * t25) * r_i_i_C(1) + (-t11 * t27 - t22 * t14) * r_i_i_C(2), 0;];
Ja_transl  = t7;
