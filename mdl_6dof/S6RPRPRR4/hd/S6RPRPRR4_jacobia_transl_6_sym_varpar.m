% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:49
% EndTime: 2019-02-26 20:50:49
% DurationCPUTime: 0.11s
% Computational Cost: add. (163->33), mult. (146->48), div. (0->0), fcn. (161->10), ass. (0->30)
t18 = sin(qJ(3));
t20 = cos(qJ(3));
t26 = pkin(3) + r_i_i_C(3) + pkin(9) + pkin(8);
t24 = t26 * t20;
t17 = sin(qJ(5));
t25 = pkin(5) * t17 + qJ(4);
t36 = -t25 * t18 - pkin(2) - t24;
t15 = qJ(1) + pkin(10);
t11 = sin(t15);
t12 = cos(t15);
t16 = qJ(5) + qJ(6);
t13 = sin(t16);
t14 = cos(t16);
t28 = t14 * t18;
t5 = -t11 * t13 + t12 * t28;
t29 = t13 * t18;
t6 = t11 * t14 + t12 * t29;
t34 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t7 = t11 * t28 + t12 * t13;
t8 = -t11 * t29 + t12 * t14;
t33 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t19 = cos(qJ(5));
t32 = pkin(5) * t19;
t31 = r_i_i_C(1) * t14;
t30 = pkin(7) + pkin(4) + t32;
t27 = t18 * t19;
t23 = r_i_i_C(1) * t13 + r_i_i_C(2) * t14 + t25;
t22 = -t26 * t18 + t23 * t20;
t9 = t20 * t13 * r_i_i_C(2);
t1 = [-sin(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t30 * t12 + t36 * t11, 0, t22 * t12, t12 * t18 (-t11 * t17 + t12 * t27) * pkin(5) + t34, t34; cos(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t30 * t11 - t36 * t12, 0, t22 * t11, t11 * t18 (t11 * t27 + t12 * t17) * pkin(5) + t33, t33; 0, 1, t23 * t18 + t24, -t20, t9 + (-t31 - t32) * t20, -t20 * t31 + t9;];
Ja_transl  = t1;
