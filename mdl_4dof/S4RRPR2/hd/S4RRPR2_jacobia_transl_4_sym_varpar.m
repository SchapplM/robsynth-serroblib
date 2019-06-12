% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RRPR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-28 15:34
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RRPR2_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_jacobia_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR2_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_jacobia_transl_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-28 15:34:22
% EndTime: 2019-05-28 15:34:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (71->13), mult. (50->14), div. (0->0), fcn. (64->6), ass. (0->13)
t23 = pkin(2) + pkin(3);
t16 = qJ(1) + qJ(2);
t14 = sin(t16);
t15 = cos(t16);
t19 = sin(qJ(4));
t20 = cos(qJ(4));
t5 = -t14 * t19 - t15 * t20;
t6 = -t14 * t20 + t15 * t19;
t22 = -t6 * r_i_i_C(1) + t5 * r_i_i_C(2);
t21 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t18 = t14 * qJ(3) + t23 * t15 - t21;
t17 = t15 * qJ(3) - t23 * t14 - t22;
t1 = [-sin(qJ(1)) * pkin(1) + t17, t17, t14, t22; cos(qJ(1)) * pkin(1) + t18, t18, -t15, t21; 0, 0, 0, 0;];
Ja_transl  = t1;
