% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PPRR5
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 11:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PPRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PPRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:16
	% EndTime: 2019-12-29 11:59:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:21
	% EndTime: 2019-12-29 11:59:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:16
	% EndTime: 2019-12-29 11:59:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:21
	% EndTime: 2019-12-29 11:59:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(3));
	t13 = cos(qJ(3));
	t15 = qJD(3) * (r_i_i_C(1) * t12 + r_i_i_C(2) * t13);
	t1 = [0, 0, -sin(pkin(6)) * t15, 0; 0, 0, cos(pkin(6)) * t15, 0; 0, 0, (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12) * qJD(3), 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 11:59:22
	% EndTime: 2019-12-29 11:59:22
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (31->20), mult. (116->44), div. (0->0), fcn. (87->6), ass. (0->18)
	t133 = cos(qJ(3));
	t130 = sin(qJ(4));
	t132 = cos(qJ(4));
	t137 = r_i_i_C(1) * t130 + r_i_i_C(2) * t132;
	t135 = t137 * t133;
	t145 = qJD(3) * t135;
	t131 = sin(qJ(3));
	t136 = r_i_i_C(1) * t132 - r_i_i_C(2) * t130 + pkin(3);
	t142 = pkin(5) + r_i_i_C(3);
	t144 = (t136 * t131 - t142 * t133) * qJD(3);
	t141 = t130 * t131;
	t140 = t131 * t132;
	t139 = qJD(3) * t131;
	t138 = qJD(4) * t133;
	t134 = qJD(4) * t137;
	t129 = cos(pkin(6));
	t128 = sin(pkin(6));
	t1 = [0, 0, (-qJD(4) * t135 - t144) * t128, -t128 * t145 + ((-t128 * t140 - t129 * t130) * r_i_i_C(1) + (t128 * t141 - t129 * t132) * r_i_i_C(2)) * qJD(4); 0, 0, (t133 * t134 + t144) * t129, t129 * t145 + ((-t128 * t130 + t129 * t140) * r_i_i_C(1) + (-t128 * t132 - t129 * t141) * r_i_i_C(2)) * qJD(4); 0, 0, t131 * t134 + (-t142 * t131 - t136 * t133) * qJD(3), (t130 * t138 + t132 * t139) * r_i_i_C(2) + (t130 * t139 - t132 * t138) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end