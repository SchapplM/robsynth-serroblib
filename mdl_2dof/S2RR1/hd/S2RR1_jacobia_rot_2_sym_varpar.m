% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S2RR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
%
% Output:
% Ja_rot [3x2]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S2RR1_jacobia_rot_2_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_jacobia_rot_2_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_jacobia_rot_2_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:11:20
% EndTime: 2019-02-26 19:11:20
% DurationCPUTime: 0.02s
% Computational Cost: add. (6->2), mult. (13->3), div. (14->4), fcn. (24->3), ass. (0->5)
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t20 = 0.1e1 / t17 ^ 2 * t19 ^ 2;
t11 = cos(atan2(0, t17));
t1 = [0, 0; (0.1e1 + t20) / (0.1e1 + 0.1e1 / t11 ^ 2 * t20) / t11, 0; 0, 1;];
Ja_rot  = t1;
