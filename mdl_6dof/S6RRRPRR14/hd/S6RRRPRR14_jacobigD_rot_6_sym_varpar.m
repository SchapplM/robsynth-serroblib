% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR14_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:40
% EndTime: 2019-02-26 22:23:40
% DurationCPUTime: 0.06s
% Computational Cost: add. (38->16), mult. (130->39), div. (0->0), fcn. (134->8), ass. (0->25)
t177 = sin(pkin(6));
t180 = sin(qJ(2));
t197 = t177 * t180;
t181 = sin(qJ(1));
t196 = t177 * t181;
t184 = cos(qJ(1));
t195 = t177 * t184;
t194 = t180 * t181;
t193 = t180 * t184;
t183 = cos(qJ(2));
t192 = t181 * t183;
t191 = t184 * t183;
t190 = qJD(1) * t177;
t182 = cos(qJ(3));
t189 = qJD(2) * t182;
t178 = cos(pkin(6));
t188 = t178 * t191 - t194;
t187 = t178 * t192 + t193;
t186 = t178 * t193 + t192;
t185 = -t178 * t194 + t191;
t179 = sin(qJ(3));
t176 = t177 * t183 * t189 + (t178 * t182 - t179 * t197) * qJD(3);
t175 = (-t186 * t179 - t182 * t195) * qJD(3) + t188 * t189 + (t179 * t196 + t185 * t182) * qJD(1);
t174 = (-t185 * t179 + t182 * t196) * qJD(3) - t187 * t189 + (t179 * t195 - t186 * t182) * qJD(1);
t1 = [0, t184 * t190, t188 * qJD(1) + t185 * qJD(2), 0, t174, t174; 0, t181 * t190, t187 * qJD(1) + t186 * qJD(2), 0, t175, t175; 0, 0, qJD(2) * t197, 0, t176, t176;];
JgD_rot  = t1;
