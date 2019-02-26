% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:45
% EndTime: 2019-02-26 20:39:45
% DurationCPUTime: 0.08s
% Computational Cost: add. (83->16), mult. (53->16), div. (0->0), fcn. (58->8), ass. (0->13)
t13 = r_i_i_C(3) + qJ(5);
t15 = pkin(4) - r_i_i_C(2);
t7 = qJ(3) + pkin(10);
t2 = sin(t7);
t4 = cos(t7);
t17 = cos(qJ(3)) * pkin(3) + t13 * t2 + t15 * t4;
t16 = pkin(2) + t17;
t14 = r_i_i_C(1) + qJ(4) + pkin(7);
t11 = -sin(qJ(3)) * pkin(3) + t13 * t4 - t15 * t2;
t8 = qJ(1) + pkin(9);
t5 = cos(t8);
t3 = sin(t8);
t1 = [-sin(qJ(1)) * pkin(1) + t14 * t5 - t16 * t3, 0, t11 * t5, t3, t5 * t2, 0; cos(qJ(1)) * pkin(1) + t14 * t3 + t16 * t5, 0, t11 * t3, -t5, t3 * t2, 0; 0, 1, t17, 0, -t4, 0;];
Ja_transl  = t1;
