% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:08:27
% EndTime: 2019-02-26 21:08:27
% DurationCPUTime: 0.11s
% Computational Cost: add. (145->31), mult. (126->46), div. (0->0), fcn. (138->10), ass. (0->29)
t19 = cos(qJ(3));
t17 = sin(qJ(3));
t30 = r_i_i_C(3) + pkin(9) + pkin(8);
t24 = t30 * t17;
t18 = cos(qJ(4));
t9 = pkin(4) * t18 + pkin(3);
t34 = t19 * t9 + pkin(2) + t24;
t14 = qJ(1) + pkin(10);
t10 = sin(t14);
t11 = cos(t14);
t15 = qJ(4) + qJ(5);
t13 = cos(t15);
t12 = sin(t15);
t29 = t12 * t19;
t5 = t10 * t29 + t11 * t13;
t28 = t13 * t19;
t6 = -t10 * t28 + t11 * t12;
t33 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = t10 * t13 - t11 * t29;
t8 = t10 * t12 + t11 * t28;
t32 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t16 = sin(qJ(4));
t31 = pkin(4) * t16;
t27 = t16 * t19;
t25 = pkin(7) + t31;
t23 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t13;
t22 = r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + t9;
t21 = -t22 * t17 + t30 * t19;
t1 = [-sin(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t25 * t11 - t34 * t10, 0, t21 * t11 (t10 * t18 - t11 * t27) * pkin(4) + t32, t32, 0; cos(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t25 * t10 + t34 * t11, 0, t21 * t10 (-t10 * t27 - t11 * t18) * pkin(4) + t33, t33, 0; 0, 1, t22 * t19 + t24 (t23 - t31) * t17, t23 * t17, 0;];
Ja_transl  = t1;
