% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:47
% EndTime: 2019-02-26 22:31:47
% DurationCPUTime: 0.09s
% Computational Cost: add. (159->22), mult. (101->24), div. (0->0), fcn. (104->8), ass. (0->20)
t40 = r_i_i_C(3) + qJ(5);
t19 = qJ(2) + qJ(3);
t16 = qJ(4) + t19;
t13 = sin(t16);
t14 = cos(t16);
t23 = (pkin(4) - r_i_i_C(2)) * t14 + t40 * t13;
t22 = pkin(3) * cos(t19) + t23;
t39 = cos(qJ(2)) * pkin(2) + t22;
t38 = t40 * t14;
t24 = -pkin(4) * t13 - pkin(3) * sin(t19);
t37 = pkin(1) + t39;
t32 = r_i_i_C(1) + pkin(9) + pkin(8) + pkin(7);
t20 = sin(qJ(1));
t31 = t20 * t13;
t21 = cos(qJ(1));
t30 = t21 * t13;
t27 = r_i_i_C(2) * t31 + t38 * t20;
t26 = r_i_i_C(2) * t30 + t38 * t21;
t25 = -sin(qJ(2)) * pkin(2) + t24;
t1 = [-t37 * t20 + t32 * t21, t25 * t21 + t26, t24 * t21 + t26, -pkin(4) * t30 + t26, t30, 0; t32 * t20 + t37 * t21, t25 * t20 + t27, t24 * t20 + t27, -pkin(4) * t31 + t27, t31, 0; 0, t39, t22, t23, -t14, 0;];
Ja_transl  = t1;
