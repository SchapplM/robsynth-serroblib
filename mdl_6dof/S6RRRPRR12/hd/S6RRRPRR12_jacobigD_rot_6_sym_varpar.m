% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR12_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:35
% EndTime: 2019-02-26 22:22:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (38->16), mult. (130->40), div. (0->0), fcn. (134->8), ass. (0->24)
t179 = sin(pkin(6));
t184 = cos(qJ(3));
t198 = t179 * t184;
t186 = cos(qJ(1));
t197 = t179 * t186;
t182 = sin(qJ(2));
t183 = sin(qJ(1));
t196 = t182 * t183;
t195 = t182 * t186;
t185 = cos(qJ(2));
t194 = t183 * t185;
t193 = t186 * t185;
t192 = qJD(1) * t179;
t181 = sin(qJ(3));
t191 = qJD(2) * t181;
t180 = cos(pkin(6));
t190 = t180 * t193 - t196;
t189 = t180 * t194 + t195;
t188 = t180 * t195 + t194;
t187 = -t180 * t196 + t193;
t178 = t179 * t185 * t191 + (t180 * t181 + t182 * t198) * qJD(3);
t177 = (-t181 * t197 + t188 * t184) * qJD(3) + t190 * t191 + (t187 * t181 - t183 * t198) * qJD(1);
t176 = (t183 * t179 * t181 + t187 * t184) * qJD(3) - t189 * t191 + (-t188 * t181 - t184 * t197) * qJD(1);
t1 = [0, t186 * t192, t190 * qJD(1) + t187 * qJD(2), 0, t176, t176; 0, t183 * t192, t189 * qJD(1) + t188 * qJD(2), 0, t177, t177; 0, 0, t179 * qJD(2) * t182, 0, t178, t178;];
JgD_rot  = t1;
