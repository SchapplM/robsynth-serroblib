% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR4_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:19
% EndTime: 2019-02-26 20:12:19
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t134 = sin(pkin(6));
t137 = sin(qJ(3));
t147 = t134 * t137;
t138 = sin(qJ(2));
t146 = t134 * t138;
t136 = cos(pkin(6));
t145 = t136 * t138;
t140 = cos(qJ(2));
t144 = t136 * t140;
t143 = qJD(2) * t137;
t133 = sin(pkin(11));
t135 = cos(pkin(11));
t142 = t133 * t140 + t135 * t145;
t141 = -t133 * t145 + t135 * t140;
t139 = cos(qJ(3));
t1 = [0, 0, t141 * qJD(2) (t133 * t147 + t141 * t139) * qJD(3) + (-t133 * t144 - t135 * t138) * t143, 0, 0; 0, 0, t142 * qJD(2) (-t135 * t147 + t142 * t139) * qJD(3) + (-t133 * t138 + t135 * t144) * t143, 0, 0; 0, 0, qJD(2) * t146, t134 * t140 * t143 + (t136 * t137 + t139 * t146) * qJD(3), 0, 0;];
JgD_rot  = t1;
