% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR7_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobia_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_transl_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:28
% EndTime: 2019-02-26 19:57:29
% DurationCPUTime: 0.09s
% Computational Cost: add. (37->15), mult. (103->31), div. (0->0), fcn. (131->10), ass. (0->18)
t11 = cos(pkin(7));
t5 = sin(pkin(14));
t7 = sin(pkin(7));
t9 = cos(pkin(14));
t15 = (r_i_i_C(1) * t5 + r_i_i_C(2) * t9) * t11 - (r_i_i_C(3) + qJ(3)) * t7;
t8 = sin(pkin(6));
t21 = t11 * t8;
t12 = cos(pkin(6));
t13 = sin(qJ(2));
t20 = t12 * t13;
t14 = cos(qJ(2));
t19 = t12 * t14;
t16 = t9 * r_i_i_C(1) - t5 * r_i_i_C(2) + pkin(2);
t10 = cos(pkin(13));
t6 = sin(pkin(13));
t3 = -t10 * t13 - t6 * t19;
t1 = t10 * t19 - t6 * t13;
t2 = [0, t16 * t3 + t15 * (-t10 * t14 + t6 * t20) t6 * t21 - t3 * t7, 0, 0, 0; 0, t16 * t1 + t15 * (-t10 * t20 - t14 * t6) -t1 * t7 - t10 * t21, 0, 0, 0; 1 (-t15 * t13 + t16 * t14) * t8, -t8 * t14 * t7 + t12 * t11, 0, 0, 0;];
Ja_transl  = t2;
