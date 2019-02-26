% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:47
% EndTime: 2019-02-26 20:26:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (57->14), mult. (37->14), div. (0->0), fcn. (41->8), ass. (0->11)
t13 = pkin(2) + r_i_i_C(3) + qJ(5) + pkin(7);
t5 = qJ(4) + pkin(10);
t1 = sin(t5);
t3 = cos(t5);
t12 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2) - sin(qJ(4)) * pkin(4);
t11 = cos(qJ(4)) * pkin(4) + r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
t10 = qJ(3) - t12;
t6 = qJ(1) + pkin(9);
t4 = cos(t6);
t2 = sin(t6);
t7 = [-sin(qJ(1)) * pkin(1) - t13 * t2 + t10 * t4, 0, t2, t11 * t2, t4, 0; cos(qJ(1)) * pkin(1) + t13 * t4 + t10 * t2, 0, -t4, -t11 * t4, t2, 0; 0, 1, 0, t12, 0, 0;];
Ja_transl  = t7;
