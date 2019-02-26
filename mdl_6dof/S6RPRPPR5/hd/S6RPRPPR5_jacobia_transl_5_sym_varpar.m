% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:24
% EndTime: 2019-02-26 20:41:24
% DurationCPUTime: 0.08s
% Computational Cost: add. (77->16), mult. (80->19), div. (0->0), fcn. (92->7), ass. (0->14)
t5 = sin(pkin(10));
t6 = cos(pkin(10));
t13 = r_i_i_C(1) * t5 + r_i_i_C(2) * t6 + qJ(4);
t14 = pkin(3) + r_i_i_C(3) + qJ(5);
t4 = pkin(9) + qJ(3);
t2 = sin(t4);
t3 = cos(t4);
t11 = t13 * t2 + t14 * t3;
t15 = cos(pkin(9)) * pkin(2) + pkin(1) + t11;
t12 = r_i_i_C(1) * t6 - r_i_i_C(2) * t5 + pkin(4) + pkin(7) + qJ(2);
t10 = t13 * t3 - t14 * t2;
t9 = cos(qJ(1));
t8 = sin(qJ(1));
t1 = [t12 * t9 - t15 * t8, t8, t10 * t9, t9 * t2, t9 * t3, 0; t12 * t8 + t15 * t9, -t9, t10 * t8, t8 * t2, t8 * t3, 0; 0, 0, t11, -t3, t2, 0;];
Ja_transl  = t1;
