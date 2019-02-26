% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:10
% EndTime: 2019-02-26 20:39:11
% DurationCPUTime: 0.08s
% Computational Cost: add. (101->19), mult. (74->20), div. (0->0), fcn. (83->10), ass. (0->15)
t10 = cos(pkin(11));
t9 = sin(pkin(11));
t16 = r_i_i_C(1) * t10 - r_i_i_C(2) * t9 + pkin(4);
t17 = r_i_i_C(3) + qJ(5);
t7 = qJ(3) + pkin(10);
t2 = sin(t7);
t4 = cos(t7);
t19 = cos(qJ(3)) * pkin(3) + t16 * t4 + t17 * t2;
t18 = pkin(2) + t19;
t15 = r_i_i_C(1) * t9 + r_i_i_C(2) * t10 + pkin(7) + qJ(4);
t13 = -sin(qJ(3)) * pkin(3) - t16 * t2 + t17 * t4;
t8 = qJ(1) + pkin(9);
t5 = cos(t8);
t3 = sin(t8);
t1 = [-sin(qJ(1)) * pkin(1) + t15 * t5 - t18 * t3, 0, t13 * t5, t3, t5 * t2, 0; cos(qJ(1)) * pkin(1) + t15 * t3 + t18 * t5, 0, t13 * t3, -t5, t3 * t2, 0; 0, 1, t19, 0, -t4, 0;];
Ja_transl  = t1;
