% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:55
% EndTime: 2019-02-26 22:03:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (127->21), mult. (81->22), div. (0->0), fcn. (86->8), ass. (0->19)
t39 = r_i_i_C(3) + qJ(5);
t19 = qJ(2) + qJ(3);
t15 = pkin(10) + t19;
t12 = sin(t15);
t13 = cos(t15);
t22 = pkin(3) * cos(t19) + (pkin(4) - r_i_i_C(2)) * t13 + t39 * t12;
t38 = cos(qJ(2)) * pkin(2) + t22;
t37 = t39 * t13;
t23 = -pkin(4) * t12 - pkin(3) * sin(t19);
t36 = pkin(1) + t38;
t31 = r_i_i_C(1) + qJ(4) + pkin(8) + pkin(7);
t20 = sin(qJ(1));
t30 = t20 * t12;
t21 = cos(qJ(1));
t29 = t21 * t12;
t26 = r_i_i_C(2) * t30 + t37 * t20;
t25 = r_i_i_C(2) * t29 + t37 * t21;
t24 = -sin(qJ(2)) * pkin(2) + t23;
t1 = [-t36 * t20 + t31 * t21, t24 * t21 + t25, t23 * t21 + t25, t20, t29, 0; t31 * t20 + t36 * t21, t24 * t20 + t26, t23 * t20 + t26, -t21, t30, 0; 0, t38, t22, 0, -t13, 0;];
Ja_transl  = t1;
