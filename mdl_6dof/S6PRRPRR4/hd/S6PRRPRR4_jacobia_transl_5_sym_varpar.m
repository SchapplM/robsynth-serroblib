% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:52
% DurationCPUTime: 0.12s
% Computational Cost: add. (139->33), mult. (362->61), div. (0->0), fcn. (465->10), ass. (0->27)
t15 = sin(pkin(6));
t19 = sin(qJ(3));
t31 = t15 * t19;
t22 = cos(qJ(3));
t30 = t15 * t22;
t17 = cos(pkin(6));
t20 = sin(qJ(2));
t29 = t17 * t20;
t23 = cos(qJ(2));
t28 = t17 * t23;
t27 = -r_i_i_C(3) - pkin(9) + pkin(8);
t18 = sin(qJ(5));
t21 = cos(qJ(5));
t26 = t18 * r_i_i_C(1) + t21 * r_i_i_C(2) + qJ(4);
t25 = t21 * r_i_i_C(1) - t18 * r_i_i_C(2) + pkin(3) + pkin(4);
t24 = t26 * t19 + t25 * t22 + pkin(2);
t16 = cos(pkin(11));
t14 = sin(pkin(11));
t10 = t17 * t19 + t20 * t30;
t9 = -t17 * t22 + t20 * t31;
t8 = -t14 * t29 + t16 * t23;
t6 = t14 * t23 + t16 * t29;
t4 = t14 * t31 + t8 * t22;
t3 = -t14 * t30 + t8 * t19;
t2 = -t16 * t31 + t6 * t22;
t1 = t16 * t30 + t6 * t19;
t5 = [0, t27 * t8 + t24 * (-t14 * t28 - t16 * t20) -t25 * t3 + t26 * t4, t3 (-t4 * t18 + t3 * t21) * r_i_i_C(1) + (-t3 * t18 - t4 * t21) * r_i_i_C(2), 0; 0, t27 * t6 + t24 * (-t14 * t20 + t16 * t28) -t25 * t1 + t26 * t2, t1 (t1 * t21 - t2 * t18) * r_i_i_C(1) + (-t1 * t18 - t2 * t21) * r_i_i_C(2), 0; 1 (t27 * t20 + t24 * t23) * t15, t26 * t10 - t25 * t9, t9 (-t10 * t18 + t9 * t21) * r_i_i_C(1) + (-t10 * t21 - t9 * t18) * r_i_i_C(2), 0;];
Ja_transl  = t5;
