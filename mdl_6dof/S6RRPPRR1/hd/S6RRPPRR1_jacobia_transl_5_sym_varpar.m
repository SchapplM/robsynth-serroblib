% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:05
% EndTime: 2019-02-26 21:28:06
% DurationCPUTime: 0.10s
% Computational Cost: add. (106->25), mult. (116->33), div. (0->0), fcn. (137->8), ass. (0->24)
t12 = qJ(2) + pkin(10);
t10 = cos(t12);
t30 = pkin(3) + pkin(4);
t9 = sin(t12);
t32 = cos(qJ(2)) * pkin(2) + t9 * qJ(4) + t30 * t10;
t31 = pkin(1) + t32;
t18 = cos(qJ(1));
t29 = t18 * t9;
t14 = sin(qJ(5));
t28 = t10 * t14;
t26 = -pkin(8) - r_i_i_C(3) + qJ(3) + pkin(7);
t16 = sin(qJ(1));
t17 = cos(qJ(5));
t21 = -t17 * t9 + t28;
t1 = t21 * t16;
t22 = t10 * t17 + t14 * t9;
t2 = t22 * t16;
t25 = -t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
t3 = -t17 * t29 + t18 * t28;
t4 = t22 * t18;
t24 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
t23 = -t22 * r_i_i_C(1) + t21 * r_i_i_C(2);
t19 = -sin(qJ(2)) * pkin(2) + qJ(4) * t10 - t30 * t9;
t5 = [-t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t31 * t16 + t26 * t18, t19 * t18 - t24, t16, t29, t24, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t26 * t16 + t31 * t18, t19 * t16 - t25, -t18, t16 * t9, t25, 0; 0, -t23 + t32, 0, -t10, t23, 0;];
Ja_transl  = t5;
