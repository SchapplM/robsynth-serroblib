% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:07
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:07:44
% EndTime: 2019-02-22 10:07:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (38->15), mult. (88->24), div. (0->0), fcn. (141->8), ass. (0->19)
t94 = cos(qJ(1));
t93 = sin(qJ(1));
t81 = sin(qJ(6));
t82 = sin(qJ(5));
t92 = t82 * t81;
t83 = cos(qJ(6));
t91 = t82 * t83;
t84 = cos(qJ(5));
t90 = t84 * t81;
t89 = t84 * t83;
t88 = cos(pkin(9));
t87 = sin(pkin(9));
t76 = -t93 * t87 - t94 * t88;
t77 = t94 * t87 - t93 * t88;
t86 = t76 * t91 + t77 * t81;
t85 = t76 * t92 - t77 * t83;
t75 = -t76 * t81 + t77 * t91;
t74 = -t76 * t83 - t77 * t92;
t1 = [t86, 0, 0, 0, t77 * t89, t74; t75, 0, 0, 0, -t76 * t89, t85; 0, 0, 0, 0, t91, t90; -t85, 0, 0, 0, -t77 * t90, -t75; t74, 0, 0, 0, t76 * t90, t86; 0, 0, 0, 0, -t92, t89; -t76 * t84, 0, 0, 0, t77 * t82, 0; -t77 * t84, 0, 0, 0, -t76 * t82, 0; 0, 0, 0, 0, -t84, 0;];
JR_rot  = t1;
