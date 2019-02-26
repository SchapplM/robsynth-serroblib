% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:27:23
% EndTime: 2019-02-26 20:27:24
% DurationCPUTime: 0.11s
% Computational Cost: add. (120->34), mult. (176->42), div. (0->0), fcn. (228->10), ass. (0->26)
t11 = qJ(4) + pkin(10);
t10 = cos(t11);
t13 = sin(qJ(6));
t15 = cos(qJ(6));
t18 = r_i_i_C(1) * t15 - r_i_i_C(2) * t13 + pkin(5);
t28 = pkin(8) + r_i_i_C(3);
t9 = sin(t11);
t20 = t28 * t9;
t31 = -t18 * t10 - t20;
t30 = pkin(1) + pkin(2);
t27 = cos(qJ(4)) * pkin(4);
t26 = cos(qJ(1));
t25 = sin(qJ(1));
t24 = t10 * t13;
t23 = t10 * t15;
t22 = cos(pkin(9));
t21 = sin(pkin(9));
t19 = t13 * r_i_i_C(1) + t15 * r_i_i_C(2);
t17 = sin(qJ(4)) * pkin(4) - t28 * t10 + t18 * t9;
t12 = -qJ(5) - pkin(7);
t8 = pkin(3) + t27;
t4 = t26 * t21 - t25 * t22;
t3 = -t25 * t21 - t26 * t22;
t2 = t13 * t4 - t3 * t23;
t1 = t15 * t4 + t3 * t24;
t5 = [t26 * qJ(2) + (-t12 + t19) * t3 + (t8 - t31) * t4 - t30 * t25, t25, 0, t17 * t3, t4, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t25 * qJ(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t4 * t12 + (-t10 * pkin(5) - t20 - t8) * t3 + t30 * t26, -t26, 0, t17 * t4, -t3 (-t15 * t3 + t4 * t24) * r_i_i_C(1) + (t13 * t3 + t4 * t23) * r_i_i_C(2); 0, 0, -1, -t27 + t31, 0, t19 * t9;];
Ja_transl  = t5;
