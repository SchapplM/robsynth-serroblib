% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(6)) * t14, 0, 0; 0, sin(pkin(6)) * t14, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->6), mult. (42->14), div. (0->0), fcn. (30->4), ass. (0->9)
	t97 = -pkin(2) + r_i_i_C(2);
	t96 = r_i_i_C(3) + qJ(3);
	t93 = cos(qJ(2));
	t95 = qJD(2) * t93;
	t92 = sin(qJ(2));
	t94 = qJD(3) * t93 + (-t96 * t92 + t97 * t93) * qJD(2);
	t91 = cos(pkin(6));
	t90 = sin(pkin(6));
	t1 = [0, t94 * t91, t91 * t95, 0; 0, t94 * t90, t90 * t95, 0; 0, t92 * qJD(3) + (t97 * t92 + t96 * t93) * qJD(2), qJD(2) * t92, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:02
	% EndTime: 2019-12-31 16:26:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (37->20), mult. (134->45), div. (0->0), fcn. (102->6), ass. (0->18)
	t131 = sin(qJ(4));
	t132 = sin(qJ(2));
	t145 = t131 * t132;
	t133 = cos(qJ(4));
	t144 = t132 * t133;
	t143 = qJD(2) * t132;
	t134 = cos(qJ(2));
	t142 = qJD(2) * t134;
	t141 = qJD(4) * t134;
	t140 = -pkin(2) - pkin(5) - r_i_i_C(3);
	t139 = r_i_i_C(1) * t133 - r_i_i_C(2) * t131;
	t138 = r_i_i_C(1) * t131 + r_i_i_C(2) * t133 + qJ(3);
	t137 = t139 * t142;
	t136 = t139 * qJD(4) + qJD(3);
	t135 = t136 * t134 + (-t138 * t132 + t140 * t134) * qJD(2);
	t130 = cos(pkin(6));
	t129 = sin(pkin(6));
	t1 = [0, t135 * t130, t130 * t142, t130 * t137 + ((-t129 * t133 - t130 * t145) * r_i_i_C(1) + (t129 * t131 - t130 * t144) * r_i_i_C(2)) * qJD(4); 0, t135 * t129, t129 * t142, t129 * t137 + ((-t129 * t145 + t130 * t133) * r_i_i_C(1) + (-t129 * t144 - t130 * t131) * r_i_i_C(2)) * qJD(4); 0, t136 * t132 + (t140 * t132 + t138 * t134) * qJD(2), t143, (-t131 * t143 + t133 * t141) * r_i_i_C(2) + (t131 * t141 + t133 * t143) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end