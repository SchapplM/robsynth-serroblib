% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:12
	% EndTime: 2019-12-31 19:53:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:12
	% EndTime: 2019-12-31 19:53:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:12
	% EndTime: 2019-12-31 19:53:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:12
	% EndTime: 2019-12-31 19:53:12
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (58->11), mult. (42->13), div. (0->0), fcn. (24->4), ass. (0->11)
	t32 = r_i_i_C(3) + qJ(3);
	t26 = qJ(1) + qJ(2);
	t23 = sin(t26);
	t25 = qJD(1) + qJD(2);
	t31 = t25 * t23;
	t24 = cos(t26);
	t30 = t25 * t24;
	t29 = pkin(1) * qJD(1);
	t28 = t23 * qJD(3) + (-pkin(2) + r_i_i_C(2)) * t31 + t32 * t30;
	t27 = r_i_i_C(2) * t30 + t24 * qJD(3) + (-pkin(2) * t24 - t32 * t23) * t25;
	t1 = [-cos(qJ(1)) * t29 + t27, t27, t30, 0, 0; -sin(qJ(1)) * t29 + t28, t28, t31, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:12
	% EndTime: 2019-12-31 19:53:12
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (117->23), mult. (118->38), div. (0->0), fcn. (74->6), ass. (0->20)
	t45 = qJ(1) + qJ(2);
	t42 = sin(t45);
	t43 = cos(t45);
	t47 = cos(qJ(4));
	t56 = qJD(4) * t47;
	t44 = qJD(1) + qJD(2);
	t46 = sin(qJ(4));
	t60 = t44 * t46;
	t62 = t42 * t56 + t43 * t60;
	t61 = t44 * t43;
	t59 = t44 * t47;
	t58 = pkin(1) * qJD(1);
	t57 = qJD(4) * t46;
	t55 = -pkin(2) - pkin(7) - r_i_i_C(3);
	t53 = t43 * t59;
	t51 = t43 * t57;
	t50 = t43 * t56;
	t49 = r_i_i_C(2) * t53 + qJ(3) * t61 + (-r_i_i_C(2) * t57 + t44 * t55 + qJD(3)) * t42 + t62 * r_i_i_C(1);
	t48 = -r_i_i_C(2) * t51 + r_i_i_C(1) * t50 + t43 * qJD(3) + (t55 * t43 + (-r_i_i_C(1) * t46 - r_i_i_C(2) * t47 - qJ(3)) * t42) * t44;
	t1 = [-cos(qJ(1)) * t58 + t48, t48, t61, -t62 * r_i_i_C(2) + (-t42 * t57 + t53) * r_i_i_C(1), 0; -sin(qJ(1)) * t58 + t49, t49, t44 * t42, (-t42 * t60 + t50) * r_i_i_C(2) + (t42 * t59 + t51) * r_i_i_C(1), 0; 0, 0, 0, (-r_i_i_C(1) * t47 + r_i_i_C(2) * t46) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:12
	% EndTime: 2019-12-31 19:53:12
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (200->28), mult. (216->36), div. (0->0), fcn. (146->6), ass. (0->26)
	t206 = pkin(4) + r_i_i_C(1);
	t205 = r_i_i_C(3) + qJ(5);
	t184 = sin(qJ(4));
	t185 = cos(qJ(4));
	t189 = t205 * t185;
	t196 = t206 * t184;
	t213 = (-t196 + t189) * qJD(4) + qJD(5) * t184;
	t183 = qJ(1) + qJ(2);
	t180 = sin(t183);
	t182 = qJD(1) + qJD(2);
	t212 = t182 * t180;
	t181 = cos(t183);
	t203 = t182 * t181;
	t201 = qJD(4) * t184;
	t202 = t182 * t185;
	t211 = t180 * t202 + t181 * t201;
	t198 = t185 * qJD(5);
	t208 = (-pkin(2) - pkin(7) - r_i_i_C(2)) * t182 + qJD(3) - t198;
	t207 = t205 * t184 + t206 * t185;
	t204 = pkin(1) * qJD(1);
	t200 = qJD(4) * t185;
	t191 = t180 * t201;
	t188 = t207 * t182;
	t187 = t208 * t180 + t205 * t191 + t206 * (t180 * t200 + t184 * t203) + (qJ(3) - t189) * t203;
	t186 = (-qJ(3) - t196) * t212 + t205 * t211 + (t206 * t200 + t208) * t181;
	t1 = [-cos(qJ(1)) * t204 + t186, t186, t203, t213 * t180 + t181 * t188, -t181 * t202 + t191; -sin(qJ(1)) * t204 + t187, t187, t212, t180 * t188 - t213 * t181, -t211; 0, 0, 0, -t207 * qJD(4) + t198, t200;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end