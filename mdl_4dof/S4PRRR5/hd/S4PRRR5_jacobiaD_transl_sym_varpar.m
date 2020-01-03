% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRRR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(7)) * t14, 0, 0; 0, sin(pkin(7)) * t14, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (33->9), mult. (40->14), div. (0->0), fcn. (25->6), ass. (0->14)
	t35 = pkin(2) * qJD(2);
	t27 = qJ(2) + qJ(3);
	t25 = cos(t27);
	t26 = qJD(2) + qJD(3);
	t34 = r_i_i_C(1) * t25 * t26;
	t24 = sin(t27);
	t33 = r_i_i_C(2) * t24 * t26;
	t32 = (-r_i_i_C(1) * t24 - r_i_i_C(2) * t25) * t26;
	t31 = -cos(qJ(2)) * t35 - t34;
	t29 = cos(pkin(7));
	t28 = sin(pkin(7));
	t23 = t29 * t33;
	t22 = t28 * t33;
	t1 = [0, t31 * t29 + t23, -t29 * t34 + t23, 0; 0, t31 * t28 + t22, -t28 * t34 + t22, 0; 0, -sin(qJ(2)) * t35 + t32, t32, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:56
	% EndTime: 2019-12-31 16:33:56
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (142->29), mult. (192->54), div. (0->0), fcn. (139->8), ass. (0->29)
	t189 = qJ(2) + qJ(3);
	t186 = sin(t189);
	t187 = cos(t189);
	t193 = cos(qJ(4));
	t207 = qJD(4) * t193;
	t188 = qJD(2) + qJD(3);
	t192 = sin(qJ(4));
	t214 = t188 * t192;
	t220 = t186 * t207 + t187 * t214;
	t219 = pkin(6) + r_i_i_C(3);
	t208 = qJD(4) * t192;
	t203 = t186 * t208;
	t218 = r_i_i_C(1) * t203 + t220 * r_i_i_C(2);
	t217 = pkin(2) * qJD(2);
	t216 = t186 * t188;
	t213 = t188 * t193;
	t190 = sin(pkin(7));
	t212 = t190 * t192;
	t211 = t190 * t193;
	t191 = cos(pkin(7));
	t210 = t191 * t192;
	t209 = t191 * t193;
	t205 = t218 * t190;
	t204 = t218 * t191;
	t198 = (r_i_i_C(1) * t192 + r_i_i_C(2) * t193) * t216;
	t197 = ((-r_i_i_C(1) * t193 - pkin(3)) * t187 - t219 * t186) * t188;
	t196 = -cos(qJ(2)) * t217 + t197;
	t195 = -pkin(3) * t216 + (-t186 * t213 - t187 * t208) * r_i_i_C(1) + t219 * t187 * t188 + (t186 * t214 - t187 * t207) * r_i_i_C(2);
	t1 = [0, t196 * t191 + t204, t191 * t197 + t204, t191 * t198 + ((-t187 * t209 - t212) * r_i_i_C(1) + (t187 * t210 - t211) * r_i_i_C(2)) * qJD(4); 0, t196 * t190 + t205, t190 * t197 + t205, t190 * t198 + ((-t187 * t211 + t210) * r_i_i_C(1) + (t187 * t212 + t209) * r_i_i_C(2)) * qJD(4); 0, -sin(qJ(2)) * t217 + t195, t195, (-t187 * t213 + t203) * r_i_i_C(2) - t220 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end