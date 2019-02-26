% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:45
% EndTime: 2019-02-26 19:48:45
% DurationCPUTime: 0.10s
% Computational Cost: add. (117->24), mult. (179->42), div. (0->0), fcn. (226->9), ass. (0->26)
t15 = pkin(11) + qJ(4);
t13 = sin(t15);
t14 = cos(t15);
t24 = r_i_i_C(3) + qJ(5);
t31 = pkin(4) - r_i_i_C(2);
t32 = t24 * t13 + t31 * t14 + cos(pkin(11)) * pkin(3) + pkin(2);
t30 = r_i_i_C(1) + pkin(8) + qJ(3);
t16 = sin(pkin(10));
t17 = sin(pkin(6));
t29 = t16 * t17;
t18 = cos(pkin(10));
t28 = t17 * t18;
t21 = sin(qJ(2));
t27 = t17 * t21;
t19 = cos(pkin(6));
t26 = t19 * t21;
t22 = cos(qJ(2));
t25 = t19 * t22;
t10 = -t16 * t26 + t18 * t22;
t9 = t16 * t25 + t18 * t21;
t8 = t16 * t22 + t18 * t26;
t7 = t16 * t21 - t18 * t25;
t5 = t13 * t27 - t14 * t19;
t3 = t10 * t13 - t14 * t29;
t1 = t13 * t8 + t14 * t28;
t2 = [0, t30 * t10 - t32 * t9, t9, t24 * (t10 * t14 + t13 * t29) - t31 * t3, t3, 0; 0, t30 * t8 - t32 * t7, t7, t24 * (-t13 * t28 + t14 * t8) - t31 * t1, t1, 0; 1 (t30 * t21 + t32 * t22) * t17, -t17 * t22, t24 * (t13 * t19 + t14 * t27) - t31 * t5, t5, 0;];
Ja_transl  = t2;
