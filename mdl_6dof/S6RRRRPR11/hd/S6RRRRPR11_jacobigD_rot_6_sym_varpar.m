% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR11_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:21
% EndTime: 2019-02-26 22:36:21
% DurationCPUTime: 0.06s
% Computational Cost: add. (38->16), mult. (130->40), div. (0->0), fcn. (134->8), ass. (0->24)
t180 = sin(pkin(6));
t185 = cos(qJ(3));
t199 = t180 * t185;
t187 = cos(qJ(1));
t198 = t180 * t187;
t183 = sin(qJ(2));
t184 = sin(qJ(1));
t197 = t183 * t184;
t196 = t183 * t187;
t186 = cos(qJ(2));
t195 = t184 * t186;
t194 = t187 * t186;
t193 = qJD(1) * t180;
t182 = sin(qJ(3));
t192 = qJD(2) * t182;
t181 = cos(pkin(6));
t191 = t181 * t194 - t197;
t190 = t181 * t195 + t196;
t189 = t181 * t196 + t195;
t188 = -t181 * t197 + t194;
t179 = t180 * t186 * t192 + (t181 * t182 + t183 * t199) * qJD(3);
t178 = (-t182 * t198 + t189 * t185) * qJD(3) + t191 * t192 + (t188 * t182 - t184 * t199) * qJD(1);
t177 = (t184 * t180 * t182 + t188 * t185) * qJD(3) - t190 * t192 + (-t189 * t182 - t185 * t198) * qJD(1);
t1 = [0, t187 * t193, t191 * qJD(1) + t188 * qJD(2), t177, 0, t177; 0, t184 * t193, t190 * qJD(1) + t189 * qJD(2), t178, 0, t178; 0, 0, t180 * qJD(2) * t183, t179, 0, t179;];
JgD_rot  = t1;
