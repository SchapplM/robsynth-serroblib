% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:21
% EndTime: 2019-02-26 21:25:21
% DurationCPUTime: 0.11s
% Computational Cost: add. (118->30), mult. (133->40), div. (0->0), fcn. (149->8), ass. (0->23)
t22 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(8);
t10 = qJ(2) + pkin(9);
t8 = cos(t10);
t31 = -t22 * t8 - cos(qJ(2)) * pkin(2);
t13 = sin(qJ(5));
t7 = sin(t10);
t30 = -(pkin(5) * t13 + qJ(4)) * t7 - pkin(1) + t31;
t28 = pkin(5) + r_i_i_C(1);
t16 = cos(qJ(5));
t27 = pkin(5) * t16 + pkin(4) + pkin(7) + qJ(3);
t17 = cos(qJ(1));
t26 = t13 * t17;
t15 = sin(qJ(1));
t25 = t15 * t13;
t24 = t15 * t16;
t23 = t16 * t17;
t1 = t7 * t23 - t25;
t3 = t7 * t24 + t26;
t19 = r_i_i_C(2) * t16 + t28 * t13 + qJ(4);
t18 = -sin(qJ(2)) * pkin(2) + t19 * t8 - t22 * t7;
t4 = -t7 * t25 + t23;
t2 = t7 * t26 + t24;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t30 * t15 + t27 * t17, t18 * t17, t15, t17 * t7, -t2 * r_i_i_C(2) + t28 * t1, t17 * t8; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t27 * t15 - t30 * t17, t18 * t15, -t17, t15 * t7, t4 * r_i_i_C(2) + t28 * t3, t15 * t8; 0, t19 * t7 - t31, 0, -t8 (r_i_i_C(2) * t13 - t28 * t16) * t8, t7;];
Ja_transl  = t5;
