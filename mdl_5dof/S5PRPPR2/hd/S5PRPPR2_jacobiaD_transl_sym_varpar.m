% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPPR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:25:09
	% EndTime: 2019-12-05 15:25:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:25:09
	% EndTime: 2019-12-05 15:25:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:25:08
	% EndTime: 2019-12-05 15:25:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(7)) * t14, 0, 0, 0; 0, sin(pkin(7)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:25:08
	% EndTime: 2019-12-05 15:25:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->5), mult. (24->10), div. (0->0), fcn. (15->6), ass. (0->5)
	t15 = qJ(2) + pkin(8);
	t13 = sin(t15);
	t14 = cos(t15);
	t19 = qJD(2) * (-cos(qJ(2)) * pkin(2) - r_i_i_C(1) * t14 + r_i_i_C(2) * t13);
	t1 = [0, cos(pkin(7)) * t19, 0, 0, 0; 0, sin(pkin(7)) * t19, 0, 0, 0; 0, (-sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:25:09
	% EndTime: 2019-12-05 15:25:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (39->10), mult. (64->18), div. (0->0), fcn. (46->8), ass. (0->10)
	t130 = r_i_i_C(3) + qJ(4);
	t121 = qJ(2) + pkin(8);
	t120 = cos(t121);
	t129 = qJD(2) * t120;
	t128 = -r_i_i_C(1) * cos(pkin(9)) + r_i_i_C(2) * sin(pkin(9)) - pkin(3);
	t119 = sin(t121);
	t127 = qJD(4) * t120 + (-cos(qJ(2)) * pkin(2) - t130 * t119 + t128 * t120) * qJD(2);
	t125 = cos(pkin(7));
	t123 = sin(pkin(7));
	t1 = [0, t127 * t125, 0, t125 * t129, 0; 0, t127 * t123, 0, t123 * t129, 0; 0, t119 * qJD(4) + (-sin(qJ(2)) * pkin(2) + t130 * t120 + t128 * t119) * qJD(2), 0, qJD(2) * t119, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:25:09
	% EndTime: 2019-12-05 15:25:09
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (110->26), mult. (137->51), div. (0->0), fcn. (105->9), ass. (0->19)
	t171 = r_i_i_C(3) + pkin(6) + qJ(4);
	t157 = qJ(2) + pkin(8);
	t155 = cos(t157);
	t158 = sin(pkin(7));
	t170 = t155 * t158;
	t159 = cos(pkin(7));
	t169 = t155 * t159;
	t168 = qJD(2) * t155;
	t153 = sin(t157);
	t167 = qJD(5) * t153;
	t156 = pkin(9) + qJ(5);
	t152 = sin(t156);
	t154 = cos(t156);
	t166 = r_i_i_C(1) * t152 + r_i_i_C(2) * t154;
	t165 = -r_i_i_C(1) * t154 + r_i_i_C(2) * t152 - cos(pkin(9)) * pkin(4) - pkin(3);
	t164 = t166 * t153;
	t163 = qJD(2) * t164;
	t162 = qJD(4) * t155 + qJD(5) * t164 + (-cos(qJ(2)) * pkin(2) - t153 * t171 + t165 * t155) * qJD(2);
	t1 = [0, t162 * t159, 0, t159 * t168, t159 * t163 + ((-t152 * t158 - t154 * t169) * r_i_i_C(1) + (t152 * t169 - t154 * t158) * r_i_i_C(2)) * qJD(5); 0, t162 * t158, 0, t158 * t168, t158 * t163 + ((t152 * t159 - t154 * t170) * r_i_i_C(1) + (t152 * t170 + t154 * t159) * r_i_i_C(2)) * qJD(5); 0, t153 * qJD(4) - t166 * t155 * qJD(5) + (-sin(qJ(2)) * pkin(2) + t171 * t155 + t165 * t153) * qJD(2), 0, qJD(2) * t153, (t152 * t167 - t154 * t168) * r_i_i_C(2) + (-t152 * t168 - t154 * t167) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end