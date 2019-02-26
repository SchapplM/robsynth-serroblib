% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:04:32
% EndTime: 2019-02-26 22:04:32
% DurationCPUTime: 0.10s
% Computational Cost: add. (92->19), mult. (87->20), div. (0->0), fcn. (92->6), ass. (0->19)
t37 = r_i_i_C(1) + qJ(4);
t15 = qJ(2) + qJ(3);
t12 = sin(t15);
t13 = cos(t15);
t32 = -pkin(3) - pkin(4);
t21 = t37 * t12 + (-r_i_i_C(2) - t32) * t13;
t36 = cos(qJ(2)) * pkin(2) + t21;
t35 = t37 * t13;
t33 = pkin(1) + t36;
t17 = sin(qJ(1));
t30 = t17 * t12;
t18 = cos(qJ(1));
t29 = t18 * t12;
t25 = r_i_i_C(2) * t30 + t35 * t17;
t24 = r_i_i_C(2) * t29 + t35 * t18;
t23 = -r_i_i_C(3) - qJ(5) + pkin(8) + pkin(7);
t22 = t32 * t12;
t20 = -sin(qJ(2)) * pkin(2) + t22;
t1 = [-t33 * t17 + t23 * t18, t20 * t18 + t24, t18 * t22 + t24, t29, -t17, 0; t23 * t17 + t33 * t18, t20 * t17 + t25, t17 * t22 + t25, t30, t18, 0; 0, t36, t21, -t13, 0, 0;];
Ja_transl  = t1;
