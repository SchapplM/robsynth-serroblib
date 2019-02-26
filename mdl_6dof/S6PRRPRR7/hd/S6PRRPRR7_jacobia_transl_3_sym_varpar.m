% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR7_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobia_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR7_jacobia_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobia_transl_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:36
% EndTime: 2019-02-26 20:07:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (35->19), mult. (92->38), div. (0->0), fcn. (112->8), ass. (0->18)
t6 = sin(pkin(6));
t9 = sin(qJ(3));
t19 = t6 * t9;
t18 = pkin(8) + r_i_i_C(3);
t11 = cos(qJ(3));
t17 = t11 * t6;
t12 = cos(qJ(2));
t8 = cos(pkin(6));
t16 = t12 * t8;
t10 = sin(qJ(2));
t5 = sin(pkin(11));
t15 = t5 * t10;
t7 = cos(pkin(11));
t14 = t7 * t10;
t13 = r_i_i_C(1) * t11 - r_i_i_C(2) * t9 + pkin(2);
t4 = t12 * t7 - t8 * t15;
t2 = t12 * t5 + t8 * t14;
t1 = [0, t18 * t4 + t13 * (-t5 * t16 - t14) (t5 * t17 - t4 * t9) * r_i_i_C(1) + (-t11 * t4 - t5 * t19) * r_i_i_C(2), 0, 0, 0; 0, t18 * t2 + t13 * (t7 * t16 - t15) (-t7 * t17 - t2 * t9) * r_i_i_C(1) + (-t11 * t2 + t7 * t19) * r_i_i_C(2), 0, 0, 0; 1 (t18 * t10 + t13 * t12) * t6 (-t10 * t19 + t11 * t8) * r_i_i_C(1) + (-t10 * t17 - t8 * t9) * r_i_i_C(2), 0, 0, 0;];
Ja_transl  = t1;
