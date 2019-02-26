% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR11_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR11_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:19
% EndTime: 2019-02-26 21:34:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (119->34), mult. (211->55), div. (0->0), fcn. (263->10), ass. (0->26)
t30 = pkin(8) + cos(pkin(11)) * pkin(4) + pkin(3);
t15 = sin(pkin(6));
t18 = sin(qJ(1));
t29 = t15 * t18;
t19 = cos(qJ(2));
t28 = t15 * t19;
t20 = cos(qJ(1));
t27 = t15 * t20;
t26 = cos(pkin(6));
t25 = -r_i_i_C(3) - pkin(9) - qJ(4) - pkin(2);
t24 = sin(pkin(11)) * pkin(4) + qJ(3);
t23 = t18 * t26;
t22 = t20 * t26;
t13 = pkin(11) + qJ(5);
t11 = sin(t13);
t12 = cos(t13);
t21 = t11 * r_i_i_C(1) + t12 * r_i_i_C(2) + t24;
t17 = sin(qJ(2));
t7 = t12 * t27;
t6 = -t17 * t23 + t20 * t19;
t5 = t20 * t17 + t19 * t23;
t4 = t17 * t22 + t18 * t19;
t3 = t18 * t17 - t19 * t22;
t2 = t5 * t11 + t12 * t29;
t1 = -t11 * t29 + t5 * t12;
t8 = [-t18 * pkin(1) + t7 * r_i_i_C(1) - t21 * t3 + (-t11 * r_i_i_C(2) + t30) * t27 + t25 * t4, t21 * t6 + t25 * t5, t5, t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t20 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t24 * t5 - t25 * t6 + t30 * t29, t21 * t4 + t25 * t3, t3, t4 (t11 * t27 + t3 * t12) * r_i_i_C(1) + (-t11 * t3 + t7) * r_i_i_C(2), 0; 0 (t21 * t17 - t25 * t19) * t15, -t28, t15 * t17 (-t26 * t11 - t12 * t28) * r_i_i_C(1) + (t11 * t28 - t26 * t12) * r_i_i_C(2), 0;];
Ja_transl  = t8;
