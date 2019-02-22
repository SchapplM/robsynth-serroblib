% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:21
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPP3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiR_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:20:51
% EndTime: 2019-02-22 11:20:52
% DurationCPUTime: 0.04s
% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t68 = sin(qJ(2));
t69 = sin(qJ(1));
t76 = t69 * t68;
t67 = pkin(9) + qJ(4);
t65 = sin(t67);
t70 = cos(qJ(2));
t75 = t70 * t65;
t66 = cos(t67);
t74 = t70 * t66;
t71 = cos(qJ(1));
t73 = t71 * t68;
t72 = t71 * t70;
t64 = t69 * t65 + t66 * t72;
t63 = -t65 * t72 + t69 * t66;
t62 = t71 * t65 - t69 * t74;
t61 = t71 * t66 + t69 * t75;
t1 = [t62, -t66 * t73, 0, t63, 0, 0; t64, -t66 * t76, 0, -t61, 0, 0; 0, t74, 0, -t68 * t65, 0, 0; t61, t65 * t73, 0, -t64, 0, 0; t63, t65 * t76, 0, t62, 0, 0; 0, -t75, 0, -t68 * t66, 0, 0; -t76, t72, 0, 0, 0, 0; t73, t69 * t70, 0, 0, 0, 0; 0, t68, 0, 0, 0, 0;];
JR_rot  = t1;
