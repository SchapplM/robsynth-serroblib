% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR2_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR2_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:22
% EndTime: 2019-02-26 19:45:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (50->19), mult. (126->32), div. (0->0), fcn. (167->8), ass. (0->19)
t23 = pkin(3) - r_i_i_C(2);
t16 = cos(pkin(6));
t18 = cos(qJ(2));
t22 = t16 * t18;
t21 = r_i_i_C(3) + qJ(4);
t11 = sin(pkin(11));
t14 = cos(pkin(11));
t17 = sin(qJ(2));
t20 = t11 * t18 + t17 * t14;
t9 = t17 * t11 - t14 * t18;
t19 = t9 * t16;
t15 = cos(pkin(10));
t13 = sin(pkin(6));
t12 = sin(pkin(10));
t8 = t20 * t16;
t6 = t9 * t13;
t4 = t12 * t19 - t15 * t20;
t2 = -t12 * t20 - t15 * t19;
t1 = [0, t23 * t4 - t21 * (t12 * t8 + t15 * t9) + (-t12 * t22 - t15 * t17) * pkin(2), t12 * t13, -t4, 0, 0; 0, t23 * t2 - t21 * (t12 * t9 - t15 * t8) + (-t12 * t17 + t15 * t22) * pkin(2), -t15 * t13, -t2, 0, 0; 1, -t23 * t6 + (pkin(2) * t18 + t20 * t21) * t13, t16, t6, 0, 0;];
Ja_transl  = t1;
