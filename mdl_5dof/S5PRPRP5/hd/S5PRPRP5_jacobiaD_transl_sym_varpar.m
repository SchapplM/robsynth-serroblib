% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRPRP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:18
	% EndTime: 2019-12-05 15:39:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:19
	% EndTime: 2019-12-05 15:39:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:18
	% EndTime: 2019-12-05 15:39:19
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
	% StartTime: 2019-12-05 15:39:19
	% EndTime: 2019-12-05 15:39:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->7), mult. (56->16), div. (0->0), fcn. (41->6), ass. (0->9)
	t123 = r_i_i_C(3) + qJ(3);
	t119 = cos(qJ(2));
	t122 = qJD(2) * t119;
	t121 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(2);
	t118 = sin(qJ(2));
	t120 = qJD(3) * t119 + (-t123 * t118 + t121 * t119) * qJD(2);
	t117 = cos(pkin(7));
	t115 = sin(pkin(7));
	t1 = [0, t120 * t117, t117 * t122, 0, 0; 0, t120 * t115, t115 * t122, 0, 0; 0, t118 * qJD(3) + (t121 * t118 + t123 * t119) * qJD(2), qJD(2) * t118, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:19
	% EndTime: 2019-12-05 15:39:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (68->23), mult. (129->49), div. (0->0), fcn. (100->7), ass. (0->18)
	t157 = r_i_i_C(3) + pkin(6) + qJ(3);
	t143 = sin(pkin(7));
	t147 = cos(qJ(2));
	t156 = t143 * t147;
	t144 = cos(pkin(7));
	t155 = t144 * t147;
	t154 = qJD(2) * t147;
	t146 = sin(qJ(2));
	t153 = qJD(4) * t146;
	t142 = pkin(8) + qJ(4);
	t140 = sin(t142);
	t141 = cos(t142);
	t152 = r_i_i_C(1) * t140 + r_i_i_C(2) * t141;
	t151 = -r_i_i_C(1) * t141 + r_i_i_C(2) * t140 - cos(pkin(8)) * pkin(3) - pkin(2);
	t150 = t152 * t146;
	t149 = qJD(2) * t150;
	t148 = qJD(3) * t147 + qJD(4) * t150 + (-t157 * t146 + t151 * t147) * qJD(2);
	t1 = [0, t148 * t144, t144 * t154, t144 * t149 + ((-t140 * t143 - t141 * t155) * r_i_i_C(1) + (t140 * t155 - t141 * t143) * r_i_i_C(2)) * qJD(4), 0; 0, t148 * t143, t143 * t154, t143 * t149 + ((t140 * t144 - t141 * t156) * r_i_i_C(1) + (t140 * t156 + t141 * t144) * r_i_i_C(2)) * qJD(4), 0; 0, t146 * qJD(3) - t152 * t147 * qJD(4) + (t151 * t146 + t157 * t147) * qJD(2), qJD(2) * t146, (t140 * t153 - t141 * t154) * r_i_i_C(2) + (-t140 * t154 - t141 * t153) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:20
	% EndTime: 2019-12-05 15:39:20
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (153->34), mult. (259->59), div. (0->0), fcn. (212->7), ass. (0->25)
	t190 = pkin(8) + qJ(4);
	t188 = sin(t190);
	t189 = cos(t190);
	t208 = r_i_i_C(3) + qJ(5);
	t210 = pkin(4) + r_i_i_C(1);
	t212 = (t210 * t188 - t208 * t189) * qJD(4) - qJD(5) * t188;
	t209 = r_i_i_C(2) + pkin(6) + qJ(3);
	t191 = sin(pkin(7));
	t195 = cos(qJ(2));
	t207 = t191 * t195;
	t192 = cos(pkin(7));
	t206 = t192 * t195;
	t194 = sin(qJ(2));
	t205 = qJD(2) * t194;
	t204 = qJD(2) * t195;
	t203 = qJD(4) * t194;
	t201 = t191 * t205;
	t200 = t192 * t205;
	t199 = -t191 * t188 - t189 * t206;
	t198 = t192 * t188 - t189 * t207;
	t197 = -t208 * t188 - t210 * t189 - cos(pkin(8)) * pkin(3) - pkin(2);
	t196 = qJD(3) * t195 + t212 * t194 + (-t209 * t194 + t197 * t195) * qJD(2);
	t185 = t199 * qJD(4) + t188 * t200;
	t183 = t198 * qJD(4) + t188 * t201;
	t1 = [0, t196 * t192, t192 * t204, -t199 * qJD(5) - t208 * (t189 * t200 + (t188 * t206 - t189 * t191) * qJD(4)) + t210 * t185, -t185; 0, t196 * t191, t191 * t204, -t198 * qJD(5) - t208 * (t189 * t201 + (t188 * t207 + t189 * t192) * qJD(4)) + t210 * t183, -t183; 0, t194 * qJD(3) - t212 * t195 + (t197 * t194 + t209 * t195) * qJD(2), t205, (-t208 * t203 - t210 * t204) * t188 + (t208 * t204 + (-t210 * qJD(4) + qJD(5)) * t194) * t189, t188 * t204 + t189 * t203;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end