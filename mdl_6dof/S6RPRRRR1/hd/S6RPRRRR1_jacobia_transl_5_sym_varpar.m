% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:50
% EndTime: 2019-02-26 21:14:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (114->15), mult. (63->20), div. (0->0), fcn. (63->10), ass. (0->16)
t14 = qJ(3) + qJ(4);
t10 = qJ(5) + t14;
t5 = sin(t10);
t6 = cos(t10);
t20 = t6 * r_i_i_C(1) - r_i_i_C(2) * t5;
t19 = t20 + pkin(4) * cos(t14);
t24 = cos(qJ(3)) * pkin(3) + t19;
t18 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
t15 = t18 - pkin(4) * sin(t14);
t21 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7);
t17 = pkin(2) + t24;
t16 = -sin(qJ(3)) * pkin(3) + t15;
t12 = qJ(1) + pkin(11);
t8 = cos(t12);
t7 = sin(t12);
t1 = [-sin(qJ(1)) * pkin(1) + t21 * t8 - t17 * t7, 0, t16 * t8, t15 * t8, t18 * t8, 0; cos(qJ(1)) * pkin(1) + t21 * t7 + t17 * t8, 0, t16 * t7, t15 * t7, t18 * t7, 0; 0, 1, t24, t19, t20, 0;];
Ja_transl  = t1;
