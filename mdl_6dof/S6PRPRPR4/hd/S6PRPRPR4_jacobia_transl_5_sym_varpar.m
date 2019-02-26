% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:15
% EndTime: 2019-02-26 19:48:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (144->27), mult. (232->46), div. (0->0), fcn. (294->11), ass. (0->28)
t17 = pkin(11) + qJ(4);
t15 = sin(t17);
t16 = cos(t17);
t18 = sin(pkin(12));
t21 = cos(pkin(12));
t27 = r_i_i_C(1) * t21 - r_i_i_C(2) * t18 + pkin(4);
t33 = r_i_i_C(3) + qJ(5);
t36 = t33 * t15 + t27 * t16 + cos(pkin(11)) * pkin(3) + pkin(2);
t19 = sin(pkin(10));
t20 = sin(pkin(6));
t35 = t19 * t20;
t23 = sin(qJ(2));
t34 = t20 * t23;
t32 = cos(pkin(6));
t31 = cos(pkin(10));
t30 = t19 * t32;
t29 = t20 * t31;
t28 = t32 * t31;
t26 = t18 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(8) + qJ(3);
t24 = cos(qJ(2));
t10 = -t23 * t30 + t31 * t24;
t9 = t31 * t23 + t24 * t30;
t8 = t19 * t24 + t23 * t28;
t7 = t19 * t23 - t24 * t28;
t5 = t15 * t34 - t32 * t16;
t3 = t10 * t15 - t16 * t35;
t1 = t8 * t15 + t16 * t29;
t2 = [0, t26 * t10 - t36 * t9, t9, t33 * (t10 * t16 + t15 * t35) - t27 * t3, t3, 0; 0, t26 * t8 - t36 * t7, t7, t33 * (-t15 * t29 + t8 * t16) - t27 * t1, t1, 0; 1 (t26 * t23 + t36 * t24) * t20, -t20 * t24, t33 * (t32 * t15 + t16 * t34) - t27 * t5, t5, 0;];
Ja_transl  = t2;
