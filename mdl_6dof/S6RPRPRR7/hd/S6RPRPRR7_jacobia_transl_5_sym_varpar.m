% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:52:12
% EndTime: 2019-02-26 20:52:12
% DurationCPUTime: 0.08s
% Computational Cost: add. (68->18), mult. (50->18), div. (0->0), fcn. (54->8), ass. (0->16)
t10 = qJ(3) + pkin(10);
t8 = qJ(5) + t10;
t5 = sin(t8);
t6 = cos(t8);
t15 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t19 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t10) - t15;
t18 = r_i_i_C(1) * t6;
t17 = r_i_i_C(2) * t5;
t16 = pkin(1) + r_i_i_C(3) + pkin(8) + qJ(4) + pkin(7);
t14 = qJ(2) - t19;
t13 = cos(qJ(1));
t12 = sin(qJ(1));
t4 = t13 * t17;
t3 = t12 * t18;
t2 = pkin(4) * cos(t10) + cos(qJ(3)) * pkin(3);
t1 = [-t12 * t16 + t14 * t13, t12, t3 + (t2 - t17) * t12, t13, -t12 * t17 + t3, 0; t14 * t12 + t13 * t16, -t13, t4 + (-t2 - t18) * t13, t12, -t13 * t18 + t4, 0; 0, 0, t19, 0, -t15, 0;];
Ja_transl  = t1;
