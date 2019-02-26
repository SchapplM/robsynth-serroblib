% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:03
% EndTime: 2019-02-26 22:16:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (135->16), mult. (69->20), div. (0->0), fcn. (71->10), ass. (0->16)
t15 = qJ(2) + qJ(3);
t11 = pkin(11) + t15;
t10 = qJ(5) + t11;
t6 = sin(t10);
t7 = cos(t10);
t23 = t7 * r_i_i_C(1) - r_i_i_C(2) * t6;
t21 = t23 + pkin(4) * cos(t11) + pkin(3) * cos(t15);
t27 = cos(qJ(2)) * pkin(2) + t21;
t22 = -r_i_i_C(1) * t6 - r_i_i_C(2) * t7;
t18 = t22 - pkin(3) * sin(t15) - pkin(4) * sin(t11);
t24 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8) + pkin(7);
t20 = pkin(1) + t27;
t19 = -sin(qJ(2)) * pkin(2) + t18;
t17 = cos(qJ(1));
t16 = sin(qJ(1));
t1 = [-t20 * t16 + t24 * t17, t19 * t17, t18 * t17, t16, t22 * t17, 0; t24 * t16 + t20 * t17, t19 * t16, t18 * t16, -t17, t22 * t16, 0; 0, t27, t21, 0, t23, 0;];
Ja_transl  = t1;
