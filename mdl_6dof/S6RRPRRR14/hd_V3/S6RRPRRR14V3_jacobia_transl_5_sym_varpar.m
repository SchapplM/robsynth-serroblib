% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14V3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14V3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobia_transl_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:04
% DurationCPUTime: 0.12s
% Computational Cost: add. (57->32), mult. (159->62), div. (0->0), fcn. (194->8), ass. (0->30)
t10 = sin(qJ(4));
t31 = t10 * r_i_i_C(3);
t11 = sin(qJ(2));
t14 = cos(qJ(4));
t30 = t11 * t14;
t12 = sin(qJ(1));
t29 = t12 * t10;
t28 = t12 * t11;
t15 = cos(qJ(2));
t27 = t14 * t15;
t16 = cos(qJ(1));
t26 = t14 * t16;
t25 = t16 * t10;
t24 = t16 * t11;
t23 = t11 * qJ(3);
t13 = cos(qJ(5));
t9 = sin(qJ(5));
t22 = r_i_i_C(1) * t13 - r_i_i_C(2) * t9;
t4 = t12 * t27 - t25;
t21 = -t13 * t4 - t9 * t28;
t20 = t13 * t28 - t4 * t9;
t19 = t13 * t15 + t9 * t30;
t18 = -t13 * t30 + t15 * t9;
t17 = t18 * r_i_i_C(1) + t19 * r_i_i_C(2) + t15 * qJ(3) - t11 * t31;
t6 = t15 * t26 + t29;
t5 = -t12 * t14 + t15 * t25;
t3 = -t15 * t29 - t26;
t2 = t6 * t13 + t9 * t24;
t1 = t13 * t24 - t6 * t9;
t7 = [t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t3 * r_i_i_C(3) - t12 * t23, t17 * t16, t24, r_i_i_C(3) * t6 - t22 * t5, r_i_i_C(1) * t1 - t2 * r_i_i_C(2), 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * r_i_i_C(3) + t16 * t23, t17 * t12, t28, r_i_i_C(3) * t4 + t22 * t3, t20 * r_i_i_C(1) + t21 * r_i_i_C(2), 0; 0 (t11 * t9 + t13 * t27) * r_i_i_C(1) + (t11 * t13 - t9 * t27) * r_i_i_C(2) + t15 * t31 + t23, -t15 (r_i_i_C(3) * t14 - t22 * t10) * t11, -t19 * r_i_i_C(1) + t18 * r_i_i_C(2), 0;];
Ja_transl  = t7;
