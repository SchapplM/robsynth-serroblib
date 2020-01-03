% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4PRRP5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:29
	% EndTime: 2019-12-31 16:29:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:29
	% EndTime: 2019-12-31 16:29:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:29
	% EndTime: 2019-12-31 16:29:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(6));
	t15 = sin(pkin(6));
	t1 = [0, -t16 * t19, 0, 0; 0, -t15 * t19, 0, 0; 0, -t20, 0, 0; 0, t16 * t20, 0, 0; 0, t15 * t20, 0, 0; 0, -t19, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:30
	% EndTime: 2019-12-31 16:29:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (18->15), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->17)
	t154 = sin(qJ(3));
	t157 = cos(qJ(2));
	t167 = t154 * t157;
	t156 = cos(qJ(3));
	t166 = t156 * t157;
	t155 = sin(qJ(2));
	t165 = qJD(2) * t155;
	t164 = qJD(2) * t157;
	t163 = qJD(3) * t155;
	t162 = qJD(3) * t157;
	t152 = sin(pkin(6));
	t161 = t152 * t165;
	t153 = cos(pkin(6));
	t160 = t153 * t165;
	t159 = t154 * t163 - t156 * t164;
	t158 = t154 * t164 + t156 * t163;
	t1 = [0, t159 * t153, t154 * t160 + (-t152 * t154 - t153 * t166) * qJD(3), 0; 0, t159 * t152, t154 * t161 + (-t152 * t166 + t153 * t154) * qJD(3), 0; 0, -t154 * t162 - t156 * t165, -t158, 0; 0, t158 * t153, t156 * t160 + (-t152 * t156 + t153 * t167) * qJD(3), 0; 0, t158 * t152, t156 * t161 + (t152 * t167 + t153 * t156) * qJD(3), 0; 0, t154 * t165 - t156 * t162, t159, 0; 0, -t160, 0, 0; 0, -t161, 0, 0; 0, t164, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:29:30
	% EndTime: 2019-12-31 16:29:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (18->15), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->17)
	t167 = sin(qJ(3));
	t170 = cos(qJ(2));
	t180 = t167 * t170;
	t169 = cos(qJ(3));
	t179 = t169 * t170;
	t168 = sin(qJ(2));
	t178 = qJD(2) * t168;
	t177 = qJD(2) * t170;
	t176 = qJD(3) * t168;
	t175 = qJD(3) * t170;
	t165 = sin(pkin(6));
	t174 = t165 * t178;
	t166 = cos(pkin(6));
	t173 = t166 * t178;
	t172 = t167 * t176 - t169 * t177;
	t171 = t167 * t177 + t169 * t176;
	t1 = [0, t172 * t166, t167 * t173 + (-t165 * t167 - t166 * t179) * qJD(3), 0; 0, t172 * t165, t167 * t174 + (-t165 * t179 + t166 * t167) * qJD(3), 0; 0, -t167 * t175 - t169 * t178, -t171, 0; 0, t171 * t166, t169 * t173 + (-t165 * t169 + t166 * t180) * qJD(3), 0; 0, t171 * t165, t169 * t174 + (t165 * t180 + t166 * t169) * qJD(3), 0; 0, t167 * t178 - t169 * t175, t172, 0; 0, -t173, 0, 0; 0, -t174, 0, 0; 0, t177, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end