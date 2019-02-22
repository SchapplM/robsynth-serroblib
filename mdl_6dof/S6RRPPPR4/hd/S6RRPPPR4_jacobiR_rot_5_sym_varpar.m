% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:08
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:08:07
% EndTime: 2019-02-22 11:08:07
% DurationCPUTime: 0.03s
% Computational Cost: add. (10->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->11)
t63 = sin(qJ(2));
t64 = sin(qJ(1));
t70 = t64 * t63;
t65 = cos(qJ(2));
t69 = t64 * t65;
t66 = cos(qJ(1));
t68 = t66 * t63;
t67 = t66 * t65;
t62 = cos(pkin(9));
t61 = sin(pkin(9));
t1 = [-t61 * t70 + t66 * t62, t61 * t67, 0, 0, 0, 0; t61 * t68 + t64 * t62, t61 * t69, 0, 0, 0, 0; 0, t63 * t61, 0, 0, 0, 0; -t69, -t68, 0, 0, 0, 0; t67, -t70, 0, 0, 0, 0; 0, t65, 0, 0, 0, 0; t66 * t61 + t62 * t70, -t62 * t67, 0, 0, 0, 0; t64 * t61 - t62 * t68, -t62 * t69, 0, 0, 0, 0; 0, -t63 * t62, 0, 0, 0, 0;];
JR_rot  = t1;
