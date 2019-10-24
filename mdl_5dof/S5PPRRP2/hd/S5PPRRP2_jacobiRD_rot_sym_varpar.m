% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:19
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PPRRP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t24 = qJD(3) * sin(pkin(7));
	t23 = qJD(3) * cos(pkin(7));
	t20 = pkin(8) + qJ(3);
	t19 = cos(t20);
	t18 = sin(t20);
	t1 = [0, 0, -t19 * t23, 0, 0; 0, 0, -t19 * t24, 0, 0; 0, 0, -qJD(3) * t18, 0, 0; 0, 0, t18 * t23, 0, 0; 0, 0, t18 * t24, 0, 0; 0, 0, -qJD(3) * t19, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (45->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t166 = sin(pkin(7));
	t168 = sin(qJ(4));
	t182 = t166 * t168;
	t169 = cos(qJ(4));
	t181 = t166 * t169;
	t167 = cos(pkin(7));
	t180 = t167 * t168;
	t179 = t167 * t169;
	t165 = pkin(8) + qJ(3);
	t163 = sin(t165);
	t178 = qJD(3) * t163;
	t177 = qJD(3) * t168;
	t176 = qJD(3) * t169;
	t175 = qJD(4) * t168;
	t174 = qJD(4) * t169;
	t173 = t166 * t178;
	t172 = t167 * t178;
	t164 = cos(t165);
	t171 = t163 * t175 - t164 * t176;
	t170 = t163 * t174 + t164 * t177;
	t1 = [0, 0, t171 * t167, t168 * t172 + (-t164 * t179 - t182) * qJD(4), 0; 0, 0, t171 * t166, t168 * t173 + (-t164 * t181 + t180) * qJD(4), 0; 0, 0, -t163 * t176 - t164 * t175, -t170, 0; 0, 0, t170 * t167, t169 * t172 + (t164 * t180 - t181) * qJD(4), 0; 0, 0, t170 * t166, t169 * t173 + (t164 * t182 + t179) * qJD(4), 0; 0, 0, t163 * t177 - t164 * t174, t171, 0; 0, 0, -t172, 0, 0; 0, 0, -t173, 0, 0; 0, 0, qJD(3) * t164, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:42
	% EndTime: 2019-10-24 10:19:42
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (45->17), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t225 = sin(pkin(7));
	t227 = sin(qJ(4));
	t241 = t225 * t227;
	t228 = cos(qJ(4));
	t240 = t225 * t228;
	t226 = cos(pkin(7));
	t239 = t226 * t227;
	t238 = t226 * t228;
	t224 = pkin(8) + qJ(3);
	t222 = sin(t224);
	t237 = qJD(3) * t222;
	t236 = qJD(3) * t227;
	t235 = qJD(3) * t228;
	t234 = qJD(4) * t227;
	t233 = qJD(4) * t228;
	t232 = t225 * t237;
	t231 = t226 * t237;
	t223 = cos(t224);
	t230 = -t222 * t234 + t223 * t235;
	t229 = -t222 * t233 - t223 * t236;
	t1 = [0, 0, -t230 * t226, t227 * t231 + (-t223 * t238 - t241) * qJD(4), 0; 0, 0, -t230 * t225, t227 * t232 + (-t223 * t240 + t239) * qJD(4), 0; 0, 0, -t222 * t235 - t223 * t234, t229, 0; 0, 0, -t231, 0, 0; 0, 0, -t232, 0, 0; 0, 0, qJD(3) * t223, 0, 0; 0, 0, t229 * t226, -t228 * t231 + (-t223 * t239 + t240) * qJD(4), 0; 0, 0, t229 * t225, -t228 * t232 + (-t223 * t241 - t238) * qJD(4), 0; 0, 0, -t222 * t236 + t223 * t233, t230, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end